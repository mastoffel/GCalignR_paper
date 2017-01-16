### Calculating error rate in peak alignment based on substances identified by mass-spektrometry
### ============================================================================================

rm(list = ls())
library(GCalignR)
## aligned data
load("ottensmann-stoffel-hoffman/Supplementary-A/data/bbim_aligned.RData")
bbim <- bbim$aligned$RT
## MS-data
bbim_ms <- read.table("ottensmann-stoffel-hoffman/Supplementary-A/data/bbim_ms.txt",sep = "\t",header = T)
## Two decimals
bbim_ms[,2:ncol(bbim_ms)] <- round(bbim_ms[,2:ncol(bbim_ms)],digits = 2)

aligned <- as.list(bbim[,2:ncol(bbim)])
ms <- as.list(bbim_ms[,3:ncol(bbim_ms)])

## delete redundant rows which contain data on unknown substances only
indices <- list()
for (i in 1:length(aligned)) {
 indices <- append(indices,list(which(aligned[[i]] %in% ms[[i]])))   
}
rows <- sort(unique(unlist(indices)))
## update the lists
aligned <- as.list(bbim[rows,2:ncol(bbim)])
ms <- as.list(bbim_ms[rows,3:ncol(bbim_ms)])

merge <- bbim[rows,2:ncol(bbim)]
merge <- bbim[rows,2:4]
for (i in 1:length(ms)) { # All samples
    for (n in 1:length(ms[[1]])) { # All substances
      if (!is.na(ms[[i]][n])) merge[which(aligned[[i]] == ms[[i]][n]),i] <- n*1000  
    }
}
