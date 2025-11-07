#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 WORKING_DIR OUTPUT_DIR TARGETS_BED [CYTOBAND_TSV]"
  exit 1
fi

WORKING_DIR="$1"
OUTPUT_DIR="$2"
TARGETS="$3" # path to CpG islands bed file
CYTOBAND_TSV="${4:-${PWD}/cytoband_color.tsv}"  # default to cwd cytoband file if not provided

mkdir -p "$OUTPUT_DIR"

# Load modules/environment for bedtools and R
# Uncomment and adapt the following lines if you use a module system
# module load bedtools/2.30.0
# source activate r_env

echo "Filtering per-sample CpG beds for CpG islands (chrX targets)..."
shopt -s nullglob
for file in "${WORKING_DIR}"/*cpg_[1-2].bed.gz; do
    echo "  processing: $file"
    base="$(basename "$file" .bed.gz)"            
    sample="$(echo "$base" | cut -d'.' -f1)"      # modify for sample name extraction
    hap="$(echo "$base" | cut -d'_' -f2)"         # modify for haplotype extraction
    out="${OUTPUT_DIR}/${sample}_hp${hap}_chrX_filtered_cpgs.bed"
    
    zcat "$file" | bedtools intersect -wa -wb -a - -b "$TARGETS" | \
      awk -v OFS='\t' -v sample="$sample" -v hap="$hap" '{print $1, $2, $3, $4, $5, $6, $10, $11, $12, $13, $14, $15, $16, $17, $18, $22, sample, hap}' > "$out"
done
shopt -u nullglob

echo "Merging filtered sample files..."
cat "${OUTPUT_DIR}"/*_chrX_filtered_cpgs.bed > "${OUTPUT_DIR}/combined_samples_chrX_filtered_cpgs.bed"

echo "Running R: compute per-island methylation and make plot..."
Rscript - <<RSCRIPT
# R part: compute averages and plot
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

out_dir <- "${OUTPUT_DIR}"
combined_f <- file.path(out_dir, "combined_samples_chrX_filtered_cpgs.bed")
cytoband_f <- "${CYTOBAND_TSV}"
avg_out_f <- file.path(out_dir, "Combined_Samples_Avg_perIsland.tsv")
svg_out <- file.path(out_dir, "Sample_chrX_methylation.svg")

if (!file.exists(combined_f)) stop("Combined file not found: ", combined_f)
if (!file.exists(cytoband_f)) stop("Cytoband file not found: ", cytoband_f)

# read combined file (no header)
bed <- read_delim(combined_f, delim = "\t", col_names = FALSE, show_col_types = FALSE)

# Try to assign the expected column names. Adjust if the field count differs.
# Expected final table used by downstream calc/plot: columns that include 'start','stop','per_mod','Nvalid_cov','cpg','sample','Haplotype'
# The original mod_kit header used in calcMeth.R:
mod_kit_header <- c("chr","start","stop","mod","score","strand","Nvalid_cov","per_mod","Nmod","Ncanon","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall","cpg","sample","Haplotype")

if (ncol(bed) < length(mod_kit_header)) {
  stop("Unexpected column count in combined file (", ncol(bed), "). Expected at least ", length(mod_kit_header))
}

# keep only first length(mod_kit_header) columns (if there are more) then name them
bed2 <- bed %>% select(1:length(mod_kit_header))
names(bed2) <- mod_kit_header

# ensure types
bed2 <- bed2 %>%
  mutate(
    start = as.numeric(start),
    stop  = as.numeric(stop),
    per_mod = as.numeric(per_mod),
    Nvalid_cov = as.numeric(Nvalid_cov),
    sample = as.character(sample),
    Haplotype = as.character(Haplotype),
    cpg = as.character(cpg)
  )

# normalize Haplotype labels (if numeric 1/2 convert to hp1/hp2; if already hp1/hp2 leave)
bed2 <- bed2 %>%
  mutate(Haplotype = case_when(
    Haplotype %in% c("1","1.0") ~ "hp1",
    Haplotype %in% c("2","2.0") ~ "hp2",
    grepl("^hp[12]\$", Haplotype) ~ Haplotype,
    TRUE ~ Haplotype
  ))

# compute per-island means per sample/haplotype
perIsland <- bed2 %>%
  group_by(chr, cpg, Haplotype, sample) %>%
  reframe(
    start = first(start),
    stop = last(stop),
    mean_per_mod = mean(per_mod, na.rm = TRUE), # convert to percentage to match original y-axis 0-100
    std = sd(per_mod, na.rm = TRUE),
    mean_coverage = mean(Nvalid_cov, na.rm = TRUE)
  ) %>%
  ungroup()

write_tsv(perIsland, avg_out_f)

# read cytoband and filter chrX
cytoband <- read_delim(cytoband_f, delim = "\t", show_col_types = FALSE)
cytoband_x <- cytoband %>% filter(chr == "chrX")

# plotting — similar to make_plots.R
p_main <- ggplot(perIsland, aes(x = start, y = mean_per_mod, fill = Haplotype, shape = Haplotype)) +
  geom_point(color = "black", size = 3) +
  scale_fill_manual(values = c("#D95F02", "#7570B3")) +
  scale_shape_manual(values = c(21, 24)) +
  facet_grid(rows = vars(sample)) +
  xlab("CpG Island Start Position (bp)") +
  ylab("Average Fraction Methylated (%)") +
  ylim(0, 100) +
  xlim(0, max(cytoband_x\$stop, na.rm = TRUE)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.major = element_line(linewidth = 0.25, linetype = "solid", colour = "grey"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", colour = "black", size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
  )

p_rect <- ggplot(cytoband_x) +
  geom_rect(aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.1, fill = color), color = "black", linewidth = 0.1, alpha = 0.5) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 30, 10)) +
  scale_fill_identity()

combined_plot <- (p_main / p_rect) + plot_layout(ncol = 1, heights = c(9, 0.5))

svglite::svglite(file = svg_out, width = 10, height = 5)
print(combined_plot)
dev.off()

cat("R: finished. outputs:\\n")
cat("  averages:", avg_out_f, "\\n")
cat("  plot:", svg_out, "\\n")
RSCRIPT

echo "Pipeline complete. Results are in: $OUTPUT_DIR"