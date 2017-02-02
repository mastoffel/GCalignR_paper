# Visualise Simulations of the PNAS Data Set
############################################

rm(list = ls())
### load packages
library(GCalignR); library(vegan); library(ggplot2); library(cowplot); library(gtable)
### source functions
source("R/ChromaSimFunctions.R"); source("R/NMDS-Functions.R");  source("R/ggplot_dual_axis.R"); source("R/ggplot_shared_x_axis.R")

### Load factors, reorder and tweak factor levels

## factors
factors <- read.csv("data/PNAS/factors.csv",sep = ",") 
names(factors)[1] <- "Ind"
row.names(factors) <- as.character(factors$Ind)
factors <- factors[,c("colony","family","age")]
factors <- factors[order(factors$family),]
factors <- factors[order(factors$colony),]
## replace 1 & by the respective beaches
factors$colony <- plyr::revalue(as.factor(factors$colony), c("1" = "SSB", "2" = "FWB"))
factors$family <- as.factor(as.character(factors$family))

### The computation of the following code takes partly long time, hence steps were saved

## List all RData sets that belong to the simulations
# x <- list.files("data/Simulations")
# 
# ## Load the data one by one
# for (i in 1:length(x)) {
#   load(paste0("data/Simulations/",x[i]))
#   ## assign SimAlign to variables D1-D10
#   assign(paste0("D",as.character(i)),SimData[["SimAlign"]])
# }
# 
# ## Combine the data to one list
# data <- list()
# data <- append(data, D1)
# data <- append(data, D2)
# data <- append(data, D3)
# data <- append(data, D4)
# data <- append(data, D5)
# data <- append(data, D6)
# data <- append(data, D7)
# data <- append(data, D8)
# data <- append(data, D9)
# data <- append(data, D10)
# 
# save(data, file = "data/combined_simulation_data.RData")
# ## remove everything, the large list 'data' is close to 1GB
# rm(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,x,SimData,data,i) # save memory

## load data 
load(file = "data/combined_simulation_data.RData")
names(data) <- paste0("S",as.character(1:length(data)))

### Do the normalization of peak ares, log + 1 transform and run NMDS, again time consuming and steps are saved

# scent <- lapply(data, scent_extract,factors=factors) # get the scent, normalised and log+1 transformed
# save(x=scent,file = "data/scent.RData")
# start <- Sys.time()
# scent_mds <- lapply(scent, myMetaMDS,factors) # MDS using vegan::metaMDS
# stop <- Sys.time()
# stop - start # took 37 Minutes
# save(x=scent_mds,file = "data/scent_mds.RData")

## load the results of nmds, and the scent data
load(file = "data/scent_mds.RData")
load(file = "data/scent.RData")

### Calculate the PERMANOVA with adonis, saved 

# start <- Sys.time()
# scent_adonis_colony <- lapply(scent,adonis_colony) # calculates the adonis stats
# stop <- Sys.time()
# stop - start # took 2 Minutes
# save(x=scent_adonis_colony,file = "data/scent_adonis_colony.RData")

## adonis results
load(file = "data/scent_adonis_colony.RData") 


### Get the results, to create summary figures
### ADONIS & SCORED SUBSTANCES

## x = Number of not randomly manipulated peaks
## y = Adonis R2
## z = Number of scored peaks
colony_df <- data.frame(x = 1 - rep(rep(seq(0,1,0.1),each = 10),10),
                        y = unlist(lapply(scent_adonis_colony, function(x) out <- x$aov.tab$R2[1])),
                        z = unlist(lapply(data,function(x) out <- x$Logfile$Aligned[["Peaks"]]))) 
## Error rate
colony_df$x <- 1 - colony_df$x  
colony_df$x <- factor(colony_df$x)
colony_df$z <-  colony_df$z - min(colony_df$z)


### Calculate the ability to correct systematic erros
### Before this step, validate the all syst. errors are corrected if no random shifts were applied

## Error rates in order of simulated Data S1:S1100
error_rates <-  rep(rep(seq(0,1,0.1),each = 10),10)
## first ten Sample are only linearly shifted
error_rates[1:10]
## Get the linear shifts of untreated chromtograms
linear_shifts <- as.data.frame(lapply(data[c(1:10,111:120,221:230)], function(x) out <- x[["Logfile"]][["LinearShift"]][1]))
row.names(linear_shifts) <- data[[1]][["Logfile"]][["LinearShift"]][[2]]

## Calculate the standard deviation for every row, i.e. every sample
test.val <- as.vector(apply(linear_shifts, 1, sd))
## FALSE shows that all are zero, assumption is validated
any(test.val != 0)

## True shifts of the chromatograms are these
applied_shifts <- linear_shifts[,1] * -1
## Extract all applied corrections
linear_shifts <- as.data.frame(lapply(data, function(x) out <- x[["Logfile"]][["LinearShift"]][1]))
row.names(linear_shifts) <- data[[1]][["Logfile"]][["LinearShift"]][[2]]
names(linear_shifts) <- names(data)
## Row names refer to samples, column names to the data set
head(linear_shifts)[1:15]
## Calculate deviations from the correct value
## Correction should cancel out the applied shifts
## Of course the deviation can be larger than 0.05 in the end!
for (i in 1:ncol(linear_shifts)) {
  linear_shifts[,i] <- applied_shifts + linear_shifts[,i]
}
# head(linear_shifts)[1:10]

## Get the proportion of uncorrected chromatograms
xy <- function(x){
  ## proportion of peaks that are fully corrected
  ## Substract false from all and divde by the sample size
(length(x) - length(x[x > 0]))/length(x)
}
linshift_df <- data.frame(x = error_rates, z = apply(linear_shifts, 2, xy))
linshift_df$x <- factor(linshift_df$x)


# Create Plots with ggplot2
p1 <- ggplot(colony_df, aes(x, z)) + 
  geom_smooth(size = 1.5,se = T,colour = "#00bfc4",alpha = 0.1,aes(group = 1)) + 
  geom_boxplot(fill = "#00bfc4",alpha = .1,size = 0.1,weight = 1) +
  labs(y = "Additionally\n scored substances",x = "")  + theme_bw() +
  theme(axis.text.y = element_text(colour = "black",hjust = 0.5,size = 9),
        axis.text.x = element_text(colour = "white"),
        axis.title.y = element_text(colour = "black",hjust = 0.5,size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1.5),
        axis.ticks.x = element_line(colour = "white")) +  scale_x_discrete(breaks = seq(0,1,0.01), labels = seq(0,1,0.01)) 

p2 <- ggplot(colony_df, aes(x, y)) +  
  geom_smooth(size = 1.5,colour = "#f8766d",se = T,alpha = 0.1,aes(group = 1)) +
  geom_boxplot(fill = "#f8766d",alpha = 0.1,size = 0.1) +
  background_grid(minor = 'none',major = 'none') +
  labs(x = "Proportion of Peaks with random errors", y = "Adonis R²") + theme_bw() +
  theme(axis.text.y = element_text(colour = "black",hjust = 0.5,size = 9),
        axis.text.x = element_text(colour = "black",hjust = 0.5,size = 9),
        axis.title.y = element_text(colour = "black",hjust = 0.5,size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1.5)) 

plot_1 <- ggplot_shared_x_axis(plot1 = p1, plot2 = p2)
plot_1

plot_2 <- ggplot(linshift_df, aes(x, z)) + 
  geom_smooth(size = 1.5,se = T,colour = "#f8766d",alpha = 0.1,aes(group = 1)) + 
  geom_boxplot(fill = "#f8766d",alpha = .1,size = 0.1,weight = 1) +
  background_grid(minor = 'none',major = 'none') +
  labs(y = "Proportion of\n corrected systematic shifts",x = "Proportion of Peaks with random errors") + theme_bw() +
  theme(axis.text = element_text(colour = "black",hjust = 0.5,size = 9),
        axis.title = element_text(colour = "black",hjust = 0.5,size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1.5),
        axis.ticks.x = element_line(colour = "white")) +  scale_y_continuous(breaks = seq(0,1,0.1), labels = seq(0,1,0.1)) 

## export the grob and save it with the dpi and dimensions you want
ggplot2::ggsave(plot_1,
                filename = "Random_error_adonis_substances.tiff",
                path = "figures",
                width = 6,
                height = 4, units = "in",
                dpi = 300)
ggplot2::ggsave(plot_2,
                filename = "Random_error_systematic_error_correction.tiff",
                path = "figures",
                width = 6,
                height = 4, units = "in",
                dpi = 300)
