require(ggplot2)
require(grid)

depth_labeller <- function(c) {
  c$mean_coverage <- sprintf("%sx", c$mean_coverage)
  c
}

plot_sim_purity_vs_f1 <- function(error_rate) {
  base_dir = "~/incoming/somatic_mutation_rate/simulations/"
  data <- read.table(sprintf("%s/purity_vs_f1_%.2f.tsv", base_dir, error_rate), header=T)
  p <- ggplot(subset(data, mean_coverage < 300), aes(purity, f1, group=model, color=model)) + 
      geom_line(size=1) + 
      theme_classic() + 
      theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size=12), strip.text=element_text(size=12)) + 
      facet_grid(mean_coverage ~ ., labeller = (mean_coverage = depth_labeller)) + 
      scale_colour_manual(values=c("#00AFBB", "#E7B800")) + 
      ggtitle(sprintf("Simulation - 5.0 muts/mb, %d%% error", error_rate * 100))
  
  
  out <- sprintf("%s/purity_vs_f1_%.2f.pdf", base_dir, error_rate)
  ggsave(out, p, width=6, height=4)
}

d_som_unphased <- function(alt, total_depth, error_rate, purity) {
  p_som_read = (1 - error_rate) * purity / 2 + error_rate * (1 - purity/2)
  return(dbinom(alt, total_depth, p_som_read))
}

d_het_unphased <- function(alt, total_depth, error_rate, purity) {
  return(dbinom(alt, total_depth, 0.5))
}

d_som_phased <- function(alt, total_depth, error_rate, purity) {
  p_som_read = (1 - error_rate) * purity + error_rate * (1 - purity)
  return(dbinom(alt, total_depth, p_som_read))
}

d_het_phased <- function(alt, total_depth, error_rate, purity) {
  return(dbinom(alt, total_depth, 1 - error_rate))
}

plot_likelihood_df <- function(df) {
  p <- ggplot(df, aes(alt_count, p, fill=variant)) +
    geom_area(position="identity", alpha=0.5) + 
    theme_classic() + 
    theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size=12), strip.text=element_text(size=12)) + 
    scale_fill_manual(values=c("#A0B1BA", "#FF1F5B"))
  return(p)
}

# reference for color palletes: 
# https://www.molecularecologist.com/2020/04/23/simple-tools-for-mastering-color-in-scientific-figures/

plot_likelihoods <- function(purity) {
 depth = 40

 er <- 0.01

 x_unphased <- seq(0, depth) 
 df_unphased <- data.frame(model = "unphased",
                          alt_count = x_unphased, 
                          somatic = d_som_unphased(x_unphased, depth, er, purity),
                          het = d_het_unphased(x_unphased, depth, er, purity))
 
 #g_unphased <- tidyr::gather(d_unphased, "variant", "p", 2:3)
 
 x_phased <- seq(0, depth / 2)
 df_phased <- data.frame(model = "phased",
                        alt_count = x_phased, 
                        somatic = d_som_phased(x_phased, depth / 2, er, purity),
                        het = d_het_phased(x_phased, depth / 2, er, purity))
 
 df <- rbind(df_unphased, df_phased)
 
 g <- tidyr::gather(df, "variant", "p", 3:4)

 p <- ggplot(g, aes(alt_count, p, fill=variant)) +
   geom_area(position="identity", alpha=0.5) + 
   theme_classic() + 
   facet_wrap(. ~ model, scales="free") +
   theme(axis.title = element_text(size=14,face="bold"), 
        axis.text = element_text(size=12), 
        strip.text=element_text(size=12),
        strip.background.y=element_rect(color = NA),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.text = element_text(size=12)) + 
   scale_fill_manual(values=c("#A0B1BA", "#FF1F5B"))
 
 q <- ggplotGrob(p)
 lg <- linesGrob(x=unit(c(0,0),"npc"), y=unit(c(0,1),"npc"), 
                 gp=gpar(col="black", lwd=2))
 
 for (k in grep("strip-l",q$layout$name)) {
   q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
 }
 
 grid.draw(q)
 base_dir = "~/incoming/somatic_mutation_rate/simulations/"
 out <- sprintf("%s/dist_plot_purity_%.2f.pdf", base_dir, purity)
 ggsave(out, q, width=10, height=4)

}

get_df <- function(depth, er, purity, model_name) {

  if(model_name == "unphased") {
    x <- seq(0, depth) 
    df <- data.frame(model = model_name,
                     purity = purity,
                     alt_count = x, 
                     somatic = d_som_unphased(x, depth, er, purity),
                     het = d_het_unphased(x, depth, er, purity))
  
    return(df);
  }
  
  if(model_name == "phased") {
    x <- seq(0, depth / 2)
    df <- data.frame(model = model_name,
                     purity = purity,
                     alt_count = x, 
                     somatic = d_som_phased(x, depth / 2, er, purity),
                     het = d_het_phased(x, depth / 2, er, purity))
    return(df)
  }
}

purity_labeller <- function(c) {
  c$purity <- sprintf("\u03b1 = %s", c$purity)
  c
}

plot_likelihoods_v2 <- function(model_name) {
  depth = 40
  er <- 0.05
  
  df_low <- get_df(depth, er, 0.25, model_name)
  df_med <- get_df(depth, er, 0.50, model_name)
  df_hi <- get_df(depth, er, 0.75, model_name)
  
  df <- rbind(df_low, df_med, df_hi)
  print(df)
  g <- tidyr::gather(df, "class", "p", 4:5)
  
  p <- ggplot(g, aes(alt_count, p, fill=class)) +
    geom_area(position="identity", alpha=0.5) + 
    theme_classic() + 
    facet_wrap(paste("purity", purity) ~ ., scales="free") +
          theme(axis.title = element_text(size=14,face="bold"), 
          axis.text = element_text(size=12), 
          strip.text=element_text(size=12),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "bottom",
          legend.text = element_text(size=12)) + 
    scale_fill_manual(values=c("#A0B1BA", "#FF1F5B"))
  
  base_dir = "~/incoming/somatic_mutation_rate/simulations/"
  out <- sprintf("%s/dist_plot_%s.pdf", base_dir, model_name)
  ggsave(out, p, width=10, height=4)
  
}
