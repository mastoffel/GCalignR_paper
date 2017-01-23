rm(list = ls())
library(GCalignR) # load GCalignR
library(vegan)
library(ggplot2)
require(cowplot)
require(gtable)
source("R/ChromaSimFunctions.R")
source("R/NMDS-Functions.R")

# Simulate two samples containing 2 shared peaks
x <- seq(from = 3,to = 7.5,by = 1.25) # time of peak max

# create a gaussian kernel to create smooth peaks
set.seed(300)
gauss <- rnorm(99999,10,1)
gauss <- hist(gauss,breaks = 300) 
kernel <- gauss$density/100 

# Dummy areas
y <- c(4029631, 4587684, 3017028,2148778) 

# create chromatograms
kernel <- c(rep(0,45),kernel,rep(0,45)) 
df <- data.frame(x = sort(x),y = y)
a_df <- df[c(1,2),]
b_df <- df[c(1,2,4),]
c_df <- df[c(1,2,4),]

a_df <- data.frame(y = unlist(lapply(a_df$y,function(x) x * kernel)),
                   x = unlist(lapply(a_df$x,function(x) seq(x - 0.5,x + 0.499,0.002))))

b_df <- data.frame(y = unlist(lapply(b_df$y,function(x) x * kernel)) - 5000,
                   x = unlist(lapply(b_df$x,function(x) seq(x + 0.2,x + 1.199,0.002))))

c_df <- data.frame(y = unlist(lapply(c_df$y,function(x) x * kernel)) - 5000,
                   x = unlist(lapply(c_df$x,function(x) seq(x - 0.5,x + 0.499,0.002))))


# smooth curves
#------------------------------------------
a_lowpass.spline <- with(a_df,smooth.spline(x,y, spar = 0.6)) ## Control spar for amount of smoothing
a2 <- predict(a_lowpass.spline, a_df$x)
a2$y[a2$y < 0] <- 0
a_df$y <- a2$y
a_df$id <- rep("a",nrow(a_df))

b_lowpass.spline <- with(b_df,smooth.spline(x,y, spar = 0.6)) 
b2 <- predict(b_lowpass.spline, b_df$x)
b2$y[b2$y < 0] <- 0
b_df$y <- b2$y
b_df$id <- rep("b",nrow(b_df))

c_lowpass.spline <- with(c_df,smooth.spline(x,y, spar = 0.6)) ## Control spar for amount of smoothing
c2 <- predict(c_lowpass.spline, c_df$x)
c2$y[c2$y < 0] <- 0
c_df$y <- c2$y
c_df$id <- rep("c",nrow(c_df))

data <- (rbind(a_df,b_df,c_df))

# ---------------------------------------

cols <- c("Reference" = "darkred","Sample" = "darkblue","Sample_transformed" = "#7570B3")

p <- ggplot2::ggplot(data = a_df,aes(x,y)) +
  geom_rect(aes(xmin=3,xmax=3.5,ymin=0,ymax=15950), fill='pink', alpha=0.02)+
  annotate("text",x=3.25,y=16100,label="min_peak2peak", fontface = "italic",col="Black",size=3,hjust=+0.5,vjust=0)+
  geom_rect(aes(xmin=4,xmax=4.5,ymin=0,ymax=18100), fill='pink', alpha=0.02)+
  annotate("text",x=4.25,y=18250,label="max_peak2peak", fontface = "italic",col="Black",size=3,hjust=+0.5,vjust=0)+
  geom_line(col = "darkred",size=1) +  geom_area(data = a_df,fill="darkred") +
  geom_line(data=b_df,col="darkblue",size=1) + geom_area(data = b_df,fill="darkblue") +
  geom_line(data = c_df,col="#7570B3",size=1) + geom_area(data = c_df,fill="#7570B3") +
  geom_segment(size=0.7,aes(x=7.45,xend=6.75,y=4000,yend=4000),arrow = arrow(length = unit(10, "points"),type = "closed", angle = 15)) +
  annotate("text",x=7.1,y=4150,label="max_linear_shift", fontface = "italic",col="Black",size=3,hjust=+0.5,vjust=0)+
  labs(title ="",x = "Retention time", y = "") + 
  theme_bw()+theme(
    axis.title.x = element_text(size = 16,vjust=0.2),
    axis.text.x  =element_blank(), 
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()) +
scale_colour_manual(name="Legend", values = c("Reference" = "darkred", "Sample" = "darkblue", "Sample_transformed" = "#7570B3")) 

p <- p +  # geom_rect(aes(xmin=6.5,xmax=7.5,ymin=13500,ymax = 16500),fill = "grey80", color = "black")+
  annotate("text",x = 7,y=16000,label = "Reference",col="darkred",size=5,hjust=+0.5,vjust=0)+
      annotate("text",x = 7,y = 15000,label="Sample",col="darkblue",size=5,hjust=+0.5,vjust=0)+
      annotate("text",x = 7,y = 14000,label="Sample, shifted",col="#7570B3",size=5,hjust=+0.5,vjust=0)
 p
ggsave(filename = "figures/ChromSim.tiff",plot = p,width = 6, height = 6,dpi = 300,units = "in")
