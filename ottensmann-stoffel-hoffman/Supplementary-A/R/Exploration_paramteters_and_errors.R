# Determination of optimal parameters by error-rate estimation

rm(list = ls())
library(GCalignR)
library(ggplot2)
library(plotly)
source("ottensmann-stoffel-hoffman/Supplementary-A/R/optimal_params.R")
source("ottensmann-stoffel-hoffman/Supplementary-A/R/error_rate.R")

### Vary alignment parameters
############################# 
results_bfla <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/bfla.txt",
               max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
save(results_bfla,file = "~/GitHub/GCalignR_paper/results_bfla.RData")

results_bbim <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/bbim.txt",
                               max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
save(results_bbim,file = "~/GitHub/GCalignR_paper/results_bbim.RData")

results_beph <- optimal_params(data = "ottensmann-stoffel-hoffman/Supplementary-A/data/beph.txt",
                               max_diff_peak2mean = seq(from = 0.01, to = 0.05,by = 0.01),rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = seq(from = 0.01, to = 0.2, by = 0.01))
save(results_beph,file = "~/GitHub/GCalignR_paper/results_beph.RData")
#############################

### Load data
### #########
load("~/GitHub/GCalignR_paper/results_bbim.RData")
load("~/GitHub/GCalignR_paper/results_beph.RData")
load("~/GitHub/GCalignR_paper/results_bfla.RData")
### #########

### Calculate error rates
errors_bbim <- data.frame(p2p = results_bbim[[2]][["p2p"]],p2m = results_bbim[[2]][["p2m"]])
errors_bbim[["error"]] <- unlist(lapply(X = results_bbim[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/bbim_ms.txt"))*100

errors_beph <- data.frame(p2p = results_beph[[2]][["p2p"]],p2m = results_beph[[2]][["p2m"]])
errors_beph[["error"]] <- unlist(lapply(X = results_beph[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/beph_ms.txt"))*100

errors_bfla <- data.frame(p2p = results_bfla[[2]][["p2p"]],p2m = results_bfla[[2]][["p2m"]])
errors_bfla[["error"]] <- unlist(lapply(X = results_bfla[[1]],error_rate,"ottensmann-stoffel-hoffman/Supplementary-A/data/bfla_ms.txt"))*100
### #####################

### Plots
### #####################

plot_bbim <- plot_ly(data = errors_bbim,x = ~p2p, y = ~p2m, z = ~error, color = ~error, colors = c("green","blue")) %>%
add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Peak2Peak'),
                        yaxis = list(title = 'Peak2Mean'),
                        zaxis = list(title = 'Misaligned peaks [%]')))

plot_bfla <- plot_ly(data = errors_bfla,x = ~p2p, y = ~p2m, z = ~error, color = ~error, colors = c("green","blue")) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Peak2Peak'),
                        yaxis = list(title = 'Peak2Mean'),
                        zaxis = list(title = 'Misaligned peaks [%]')))

plot_beph <- plot_ly(data = errors_beph,x = ~p2p, y = ~p2m, z = ~error, color = ~error, colors = c("green","blue"),) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'Peak2Peak'),
                        yaxis = list(title = 'Peak2Mean'),
                        zaxis = list(title = 'Misaligned peaks [%]')))

### #####################
plot_bbim
plot_beph
plot_bfla
