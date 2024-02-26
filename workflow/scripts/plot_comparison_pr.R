require(ggplot2)
require(tidyr)
require(dplyr)
require(yardstick)
require(argparse)

get_pr_stats <- function(filename, vaf_threshold, label) {
  data <- read.table(filename, header=T)
  sdata <- subset(data, min_vaf >= vaf_threshold)
  p <- pr_curve(subset(sdata, in_high_confidence_region == 1), as.factor(in_truth), obs_somatic_qual, event_level="second")
  p$caller = label
  return(p)
}

make_pr_plot <- function(data, colors, sample, plot_depth, plot_regions, plot_vaf, clean) {
  # remove the infinite threshold rows
  data_clipped <- data %>% filter_at(vars(.threshold), all_vars(!is.infinite(.) & . > 0))

  gl <- element_line(color = "lightgray", linewidth = 0.25, linetype = 2)
  plot <- ggplot(data_clipped, aes(x = recall, y = precision, color=caller, linetype=strategy)) +
    geom_path(linewidth=1.5) +
    scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, by = 0.1)) +
    theme_classic() +
    theme(plot.title = element_text(size=18), 
          axis.title = element_text(size=20, face="bold"), 
          axis.text = element_text(size=12), 
          strip.text = element_text(size=12),
          legend.text = element_text(size=20),
          legend.title = element_text(size=18),
          panel.grid = gl,
          panel.grid.major = gl) +
    #xlim(0, 1) + ylim(0, 1) +
    scale_colour_manual(values=colors) + 
    scale_linetype_manual(values=c(2, 1))

  print(sample)
  if(clean) {
    plot <- plot + guides(color="none", linetype="none")
  } else {
    plot <- plot + ggtitle(sprintf("%s T: %sx, N:%sx\n%s regions\nTruth VAF >= %.0f%%", sample, plot_depth, plot_depth, plot_regions, plot_vaf * 100.0))
  }
  plot
}

if(!interactive()) {
  
  # Main
  parser <- ArgumentParser(description="Plot precision/recall curves for manuscript")
  
  parser$add_argument("--output", "-o", default="precision_recall_plot_publication_comparison.pdf",
                      help="name of the output plot")
  #parser$add_argument("--fn_pattern", "-f", required=TRUE, help="pattern of annotated call set")
  parser$add_argument("--regions", "-r", required=TRUE, help="which region set to use (GIAB vs GIAB-Called)")
  parser$add_argument("--depth", "-d", required=TRUE, help="per-sample depth")
  parser$add_argument("--sample", "-s", required=TRUE, help="sample name to plot")
  parser$add_argument("--type", "-t", required=TRUE)
  parser$add_argument("--min-vaf", "-m", required=TRUE, help="only included sites with VAF >= this value in truth set")
  parser$add_argument("--clean", "-c", required=FALSE, action="store_true", default=TRUE, help="suppress legends and title")
  
  args <- parser$parse_args()
  depth = as.numeric(args$depth)
  region_path = args$regions

  print(args$sample)
  print(region_path)

  fnp = sprintf("analysis/annotated/%s_%d_%sBL_%d_rep1.$s.annotated-%s.tsv", args$sample, depth, args$sample, depth, region_path)
  fnp = gsub("$", "%", fnp, fixed=TRUE)

  clairS_ont <- get_pr_stats(sprintf(fnp, "clairS.ont"), args$min_vaf, "clairS (ONT)")
  smrest_ont <- get_pr_stats(sprintf(fnp, "smrest.ont"), args$min_vaf, "smrest (ONT)")
  
  clairS_pacbio <- get_pr_stats(sprintf(fnp, "clairS.pacbio"), args$min_vaf, "clairS (PacBio)")
  smrest_pacbio <- get_pr_stats(sprintf(fnp, "smrest.pacbio"), args$min_vaf, "smrest (PacBio)")
  
  mutect2_to <- get_pr_stats(sprintf(fnp, "mutect2_to.illumina"), args$min_vaf, "mutect2-TO (ILMN)")
  mutect2_tn <- get_pr_stats(sprintf(fnp, "mutect2_tn.illumina"), args$min_vaf, "mutect2-TN (ILMN)")

  # add some metadata   
  clairS_ont$tech <- "ONT"
  smrest_ont$tech <- "ONT"
  
  clairS_pacbio$tech <- "Pacbio"
  smrest_pacbio$tech <- "Pacbio"
  
  mutect2_to$tech <- "ILMN"
  mutect2_tn$tech <- "ILMN"
  
  clairS_ont$strategy <- "TN"
  smrest_ont$strategy <- "TO"

  clairS_pacbio$strategy <- "TN"
  smrest_pacbio$strategy <- "TO"
  
  mutect2_to$strategy <- "TO"
  mutect2_tn$strategy <- "TN"
  
  if(args$type == "LR-vs-SR") {
      data <- rbind(clairS_ont, smrest_ont, mutect2_to, mutect2_tn)
      colors <- c("#009ADE", "#AF58BA", "#FFC61E", "#00CD6C")
  } else {
      data <- rbind(clairS_ont, smrest_ont, clairS_pacbio, smrest_pacbio)
      colors <- c("#009ADE", "#F28522", "#00CD6C", "#FF1F5B")
  }
  p <- make_pr_plot(data, colors, args$sample, depth, args$regions, as.numeric(args$min_vaf), args$clean)
  
  w <- 5
  h <- 4
  if(args$clean) {
    h <- w
  }

  ggsave(args$output, p, width=w, height=h)
}
