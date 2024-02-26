require(ggplot2)
require(tidyr)
require(argparse)

make_tmb_plot <- function(data, lo, hi, output, blood_only) {
  d = subset(data, purity == 1.0 | (purity >= lo & purity < hi))
  
  g = rev(gray.colors(10))
  cliveome_color = "red"
  ymax = 10

  if(blood_only) {
      d <- subset(d, purity == 1.0)
      colors = c(cliveome_color)
      ymax = 1.0
  } else {
      d <- subset(d, purity < 1.0)
      colors = c(g[3], g[4], g[5], g[6], g[7], cliveome_color)
      ymax = 16
  }
  gl <- element_line(color = "lightgray", linewidth = 0.25, linetype = 2)
  p <- ggplot(subset(d, total_depth > 2), aes(total_depth, hard_calls_per_mb, color=sample)) + 
    geom_point(size=2) + 
    geom_line(linewidth=1.25) + 
    theme_classic() + 
    ylab("TMB (muts/MB)") + 
    xlab("Sequencing Depth (X)") +
    xlim(0, 82) +
    ylim(0, ymax) +
    scale_colour_manual(values=colors) +
    theme_classic() +
    theme(plot.title = element_text(size=18), 
          axis.title = element_text(size=20, face="bold"), 
          axis.text = element_text(size=18), 
          strip.text = element_text(size=12),
          legend.text = element_text(size=26),
          legend.title = element_text(size=24, face="bold"),
          panel.grid = gl,
          panel.grid.major = gl)
  if(! blood_only ) {
    p <- p + geom_hline(yintercept = 14.0, linetype="dashed", color="lightblue", size=1) +
             annotate("text", x = 2, y = 14.4, label="Expected")
    p <- p + labs(color='COLO829 purity')
  }
  ggsave(output, p, width=9, height=6) 
}

if(!interactive()) {
  
  # Main
  parser <- ArgumentParser(description="Plot TMB")
  
  parser$add_argument("--output", "-o", default="tmb.pdf", help="name of the output plot")
  parser$add_argument("--tmb-file", "-f", required=TRUE, help="file containing the summarized tmb estimates")
  parser$add_argument("--blood-only", "-b", required=FALSE, action="store_true", default=FALSE, help="only plot the healthy blood sample")
  
  args <- parser$parse_args()
  data <- read.table(args$tmb_file, sep='\t', header=T)

  data$sample = sprintf("%.0f%%", 100 * data$purity)
  data$sample[data$sample == "100%"] = "Healthy Blood"
  make_tmb_plot(data, 0.3, 0.75, args$output, args$blood_only)
}
