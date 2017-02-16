## Ordiplot
library(ggplot2)
a <- ggplot(data = scent_nmds,aes(MDS1,MDS2,color = colony)) +
    geom_point(size = 2) + 
    theme_void() + 
    scale_color_manual(values = c("blue","red"), 
                      name="Breeding colony",
                      breaks=c("SSB", "FWB"),
                      labels=c("Special study beach", 
                               "Freshwater beach")) +
    theme(panel.background = element_rect(colour = "black", size = 1.25,fill = NA),aspect.ratio = 1,legend.justification=c(0,0), legend.position=c(0.01,0.01), legend.title = element_text(size=12, face="bold"),legend.background = element_rect(color = "black",fill = NA, size = 1, linetype = "solid")) 

ggplot2::ggsave(a,filename = "ottensmann-stoffel-hoffman/figures/ordiplot.png",width = 5,height = 5,units = "in",dpi = 300)

## Heatmap
library(ggplot2)
library(GCalignR)
a <- gc_heatmap(aligned_peak_data,type = "binary",threshold = 0.05,main_title = "")
a <- a + theme(plot.background = element_rect(fill = "white"))

ggplot2::ggsave(a,filename = "ottensmann-stoffel-hoffman/figures/heatmap.tiff",width = 6,height = 5,units = "in",dpi = 300)

## diagnostic plots
library(GCalignR)
# 1500 x 500
plot_1 <- plot(x = aligned_peak_data, which_plot ="peak_num",cex.lab = 1.5,cex.names = 0.8)

# 500 x 500
plot_2 <- plot(x = aligned_peak_data, which_plot ="shifts",cex.lab = 1.5,cex.names = 0.8)

plot_3 <- plot(x = aligned_peak_data, which_plot ="variation",cex.lab = 1.5,cex.names = 0.8)

plot_4 <- plot(x = aligned_peak_data, which_plot ="peaks_s",cex.lab = 1.5,cex.names = 0.8)

## Simualtion 
# 1000x1000px
with(errors_bbim,scatter3D(
    x = p2p,
    y = p2m,
    z = error, 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    xlab = "",
    ylab = "",
    zlab = "",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.axis = 1,
    zlim = c(0,0.2)))

with(errors_beph,scatter3D(
    x = p2p,
    y = p2m,
    z = error, 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    xlab = "",
    ylab = "",
    zlab = "",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.axis = 1,
    zlim = c(0,0.2)))

with(errors_bfla,scatter3D(
    x = p2p,
    y = p2m,
    z = error, 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    xlab = "",
    ylab = "",
    zlab = "",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.axis = 1,
    zlim = c(0,0.2)))
