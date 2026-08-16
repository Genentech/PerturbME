# conda environment mageck
library(MAGeCKFlute)
library(ggplot2)

figures_pn = '/ahg/regevdata/projects/PerturbCITE_ICR/202008_full_exp/figures/manuscript_figures/Fig1'
mageck_out_pn = '/ahg/regevdata/projects/PerturbCITE_ICR/202008_full_exp/guide_seq/mageck_out'

low_input_sgRNA_file = file.path(mageck_out_pn, "Low_vs_Input.sgrna_summary.txt")
low_input_sdata = read.delim(low_input_sgRNA_file, check.names = FALSE)
low_input_savepn = file.path(figures_pn, "1F_low_input_rankView.pdf")
sgRankView(low_input_sdata, top = 15, bottom = 5, filename=low_input_savepn)

high_input_sgRNA_file = file.path(mageck_out_pn, "High_vs_Input.sgrna_summary.txt")
high_input_sdata = read.delim(high_input_sgRNA_file, check.names = FALSE)
high_input_savepn = file.path(figures_pn, "1F_high_input_rankView.pdf")
sgRankView(high_input_sdata, top = 15, bottom = 5, filename=high_input_savepn)