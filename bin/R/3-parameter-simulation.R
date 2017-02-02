rm(list = ls())
library(GCalignR)


pm <- c(0.01,0.02,0.03,0.04,0.05) #peak2mean
pp <- c(0.03,0.05,0.07,0.1) #peak2peak


for (i in 1:length(pm)) { 
    for (j in 1:length(pp)) {
    assign(paste0("bb_pp_",as.character(pp),"_pm_",as.character(pm[i])),align_chromatograms(data = "data/bumblee_bee_input.txt",rt_col_name = "RT",conc_col_name = "Area",min_diff_peak2peak = pp[i],max_diff_peak2mean = pm[i]))
    }
}



