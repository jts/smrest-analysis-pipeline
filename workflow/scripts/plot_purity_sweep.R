require(ggplot2)
require(tidyr)
require(argparse)

depth_labeller <- function(c) {
  print(c)
  sprintf("%sx", c)
}

if(!interactive()) {
  
  # Main
  parser <- ArgumentParser(description="Plot accuracy metric vs mixture purity")
  
  parser$add_argument("--output", "-o", default="purity_mixture_accuracy.pdf",
                      help="name of the output plot")
  parser$add_argument("--accuracy-file", "-f", required=TRUE, help="file containing the summarized accuracy values")
  parser$add_argument("--regions", "-r", required=TRUE, help="which region set to use")
  parser$add_argument("--metric", "-m", required=TRUE, help="which accuracy measure to plot (f1, sensitivity, precision")
  parser$add_argument("--sample", "-s", required=TRUE, help="sample name to plot")
  
  args <- parser$parse_args()
  data <- read.table(args$accuracy_file, sep='\t', header=T)
  data_wide <- pivot_longer(data, cols = c("f1", "sensitivity", "precision"), names_to="metric")
  print(head(data_wide))
  subdata <- subset(data_wide, (total_depth == 40 | total_depth == 80) & technology != "pacbio" & callset_region == "annotated-HC-ont-called")
  print(subdata)
  if(args$metric == "f1") {
    subdata <- subset(subdata, metric == "f1")
  } else {
    subdata <- subset(subdata, metric != "f1")
  }
  gl <- element_line(color = "lightgray", linewidth = 0.25, linetype = 2)

  p <- ggplot(subdata, aes_string("purity", "value", color="caller")) +
    geom_line(size=0.5, linetype = 2) +
    geom_point(size = 0.75) +
    theme_classic() +
    facet_grid(total_depth ~ metric, labeller = labeller(total_depth = depth_labeller, metric = label_value)) +
    xlim(0, 1.00) + 
    ylim(0, 1.00) + 
    theme(axis.title = element_text(size=14,face="bold"), 
          axis.text = element_text(size=12), 
          strip.text = element_text(size=12),
          panel.grid = gl,
          panel.grid.major = gl) +
    scale_colour_manual(values=c("#009ADE", "#AF58BA", "#FFC61E", "#00CD6C")) +
    ggtitle(sprintf("%s - %sx depth, %s regions", args$sample, args$total_depth, args$regions))

  ggsave(args$output, p, width=8, height=4)
}
