ggplot_shared_x_axis <- function(plot1 = p1, plot2 = p2){
  ## requires ggplot
  library(ggplot2)
  ## require cowplot for switching axes
  library(cowplot)
  library(grid)
  ## swith x-axis of the upper plot to the top
  g1 <- switch_axis_position(plot1, "x")
  g2 <- ggplotGrob(plot2)
  ## set bottom row to 0 height
  g1 <- gtable_squash_rows(table = g1, rows = length(g1$height))
  ## set top two rows to 0 height
  g2 <- gtable_squash_rows(table =  g2, rows = 1:2)
  ## Combine and align the plots
  plot_grid(g1, g2, ncol = 1,align = "v")
}
