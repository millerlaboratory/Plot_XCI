# Plot_XCI
A simple script to plot XCI patterns in 46,XX samples

Using an input of compressed, haplotype-aware bedmethyl pileups generated with modkit (https://github.com/nanoporetech/modkit), this script will filter for CpGs within informative CpG islands on chrX, calculate the average methylatuon per island and then plot these values per sample.

Adjust filepaths and nameing schemes in the shell script to handle your data and use the CpG island bed file and cytoband color tsv.
