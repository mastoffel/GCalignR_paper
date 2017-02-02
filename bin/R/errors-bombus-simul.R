# # Debugger ...
# bfla_chroma$BFLA01$RT[1:4] # Original RT
# bfla_shifted$BFLA01$RT[1:4] # linear shift
# bfla_shifts[1,] # shift size
# temp$BFLA01$RT[1:4] # random errors applied
# aligned$aligned$RT$BFLA01[aligned$aligned$RT$BFLA01 > 0][1:4] # only allowed to contain linear shift by aligner
# 
# 
# aligned2 <- correct_rt(d1 = aligned,d2 = bfla_shifts2,rt_col_name = "RT")
# aligned2$aligned$RT$BFLA01[aligned2$aligned$RT$BFLA01 > 0][1:4]
# 


rm(list = ls())
#library(GCalignR)
library(ggplot2)
source("ottensmann-stoffel-hoffman/Supplementary-A/R/error_rate.R")
load("ottensmann-stoffel-hoffman/Supplementary-A/data/bfla_simulations.RData")
errors_bfla <- data.frame(noise = bfla_simulations[["noise"]])
errors_bfla[["error"]] <- unlist(lapply(X = bfla_simulations[["SimAlign"]],error_rate,Reference = "ottensmann-stoffel-hoffman/Supplementary-A/data/bfla_ms.txt",rt_col_name = "RT",linshift = FALSE))

ggplot(data = errors_bfla,aes(x = noise, y = error)) + geom_count()

dev.off()
with(errors_bfla,plot(error~as.factor(noise)))
