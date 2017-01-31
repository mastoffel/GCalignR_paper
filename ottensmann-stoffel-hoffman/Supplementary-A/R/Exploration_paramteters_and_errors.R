# Determination of optimal parameters by error-rate estimation

rm(list = ls())
library(GCalignR)
library(ggplot2)
library(plotly)
library(plot3D)
source("ottensmann-stoffel-hoffman/Supplementary-A/R/optimal_params.R")
source("ottensmann-stoffel-hoffman/Supplementary-A/R/error_rate.R")

### Vary alignment parameters
############################# 
# results_bfla <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/bfla.txt",
#                max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
# save(results_bfla,file = "~/GitHub/GCalignR_paper/results_bfla.RData")
# 
# results_bbim <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/bbim.txt",
#                                max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
# save(results_bbim,file = "~/GitHub/GCalignR_paper/results_bbim.RData")
# 
# results_beph <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/beph.txt",
#                                max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
# save(results_beph,file = "~/GitHub/GCalignR_paper/results_beph.RData")
#############################

### Load data
### #########
load("~/GitHub/GCalignR_paper/results_bbim.RData")
load("~/GitHub/GCalignR_paper/results_beph.RData")
load("~/GitHub/GCalignR_paper/results_bfla.RData")
### #########

### Calculate error rates
errors_bbim <- data.frame(p2p = results_bbim[[2]][["p2p"]],p2m = results_bbim[[2]][["p2m"]])
errors_bbim[["error"]] <- unlist(lapply(X = results_bbim[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/bbim_ms.txt"))

errors_beph <- data.frame(p2p = results_beph[[2]][["p2p"]],p2m = results_beph[[2]][["p2m"]])
errors_beph[["error"]] <- unlist(lapply(X = results_beph[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/beph_ms.txt"))

errors_bfla <- data.frame(p2p = results_bfla[[2]][["p2p"]],p2m = results_bfla[[2]][["p2m"]])
errors_bfla[["error"]] <- unlist(lapply(X = results_bfla[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/bfla_ms.txt"))
### #####################

### Remove files
rm(list = c("results_bbim", "results_beph", "results_bfla"))
### ############

### Plots
### Customise: Axis, Background-Grid...
### #####################

# Plot_ly interactive plot, nice for exploration
# plot_bbim <- plot_ly(data = errors_bbim,x = ~p2p, y = ~p2m, z = ~error, color = ~error, colors = c("green","blue")) %>%
# add_markers() %>%
#     layout(scene = list(xaxis = list(title = 'Peak2Peak'),
#                         yaxis = list(title = 'Peak2Mean'),
#                         zaxis = list(title = 'Misaligned peaks [%]')))

# plot3D
par(mfrow = c(1,3),family = "serif", mai = c(0.1,0.3,0.5,0.15))

scatter3D(
    x = errors_bbim[["p2p"]],
    y = errors_bbim[["p2m"]],
    z = errors_bbim[["error"]], 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    main = "Bombus bimaculatus\nn = 717",
    xlab = "min_diff_peak2peak",
    ylab = "max_diff_peak2mean",
    zlab = "Error rate",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.lab = 1.25,
    cex.axis = 1.25,
    cex.main = 3,
    zlim = c(0,0.3))

with(errors_beph, scatter3D(
    x = p2p,
    y = p2m,
    z = error, 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    main = "Bombus ephippiatus\nn = 782",
    xlab = "min_diff_peak2peak",
    ylab = "max_diff_peak2mean",
    zlab = "Error rate",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.lab = 1.25,
    cex.axis = 1.25,
    cex.main = 3,
    zlim = c(0,0.3)))

with(errors_bfla, scatter3D(
    x = p2p ,
    y = p2m,
    z = error, 
    pch = 19,
    size = 2,
    theta = 30,
    phi = 0,
    ticktype = "detailed",
    main = "Bombus flavifrons\nn = 457",
    xlab = "min_diff_peak2peak",
    ylab = "max_diff_peak2mean",
    zlab = "Error rate",
    bty = "g",
    colkey = FALSE,
    cex = 1.5,
    cex.lab = 1.25,
    cex.axis = 1.25,
    cex.main = 3,
    zlim = c(0,0.3))) 

### #####################

### Get best parameters
### ###################

df <- data.frame(p2p = errors_bbim[["p2p"]], p2m = errors_bbim[["p2m"]], bbim = errors_bbim[["error"]],beph = errors_beph[["error"]], bfla = errors_bfla[["error"]])
x <- function(df) mean(df[3:5])
df[["mean"]] <- apply(df,1,FUN = x)
head(df)
# Best for mean error of all three
df[which(df[["mean"]] == min(df[["mean"]])),]

# best for bbim
df[which(df[["bbim"]] == min(df[["bbim"]])),]
df[which(df[["beph"]] == min(df[["beph"]])),]
df[which(df[["bfla"]] == min(df[["bfla"]])),]