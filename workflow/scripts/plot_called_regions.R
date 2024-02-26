require(ggplot2)
require(argparse)

if(!interactive()) {
  
  # Main
  parser <- ArgumentParser(description="Plot accuracy metric vs mixture purity")
  
  parser$add_argument("--output", "-o", default="purity_mixture_accuracy.pdf",
                      help="name of the output plot")
  parser$add_argument("--called-file", "-f", required=TRUE, help="file containing the summarized called regions")
  parser$add_argument("--total-depth", "-d", required=TRUE, help="which total depth to plot")
  
  args <- parser$parse_args()
  data <- read.table(args$called_file, sep='\t', header=T)
  data$span_mb = data$span / 1000000

  subdata <- subset(data, total_depth == args$total_depth)

  p <- ggplot(subdata, aes(purity, span_mb, color=callset_region)) +
    geom_line(size=1) +
    theme_classic() +
    theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size=12), strip.text=element_text(size=12)) +
    scale_colour_manual(values=c("#00CD6C", "#AF58BA")) +
    ylab("Callset Span (Mbp)") +
    ylim(0, 2500) +
    ggtitle(sprintf("COLO829 - %sx depth", args$total_depth))

  ggsave(args$output, p, width=6, height=4)
}
