require(ggplot2)
require(tidyr)
require(dplyr)
require(yardstick)
require(argparse)

if(!interactive()) {
  
  # Main
  parser <- ArgumentParser(description="Plot precision/recall curves for manuscript")
  
  parser$add_argument("--output", "-o", default="precision_recall_plot_publication_comparison.pdf",
                      help="name of the output plot")
  parser$add_argument("--min-vaf", "-m", required=TRUE, help="annotate VAF cutoff")
  parser$add_argument("--sample", "-s", required=TRUE, help="sample name to plot")
  
  args <- parser$parse_args()
  depth = as.numeric(args$depth)
  region_path = args$regions

  fnp = sprintf("analysis/annotated/%s_40_%sBL_40_rep1.smrest.ont.annotated-HC.tsv", args$sample, args$sample)
  smrest_ont <- read.table(fnp, header=T)
  sdata <- subset(smrest_ont, in_high_confidence_region & in_truth)
  print(head(sdata))
  plot <- ggplot(sdata, aes(as.numeric(truth_vaf) * 100)) +
    geom_histogram(binwidth=1, color = "#000000", fill="darkgrey") +
    theme_classic() +
    xlab("VAF (%)") +
    geom_vline(xintercept = as.numeric(args$min_vaf) * 100, linetype="dashed", color="grey", size=1.5) +
    theme(plot.title = element_text(size=30), 
          axis.title = element_text(size=30, face="bold"), 
          axis.text = element_text(size=18), 
          strip.text = element_text(size=12),
          legend.text = element_text(size=20),
          legend.title = element_text(size=18)) + ggtitle(sprintf("%s", args$sample))

  ggsave(args$output, plot, width=10, height=8)
}
