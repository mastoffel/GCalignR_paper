### Error estimation

# The detection of errors in the alignment of gas-chromatography is not a trivial case. With the aid of additional information obtained by mass-spektrometry it is possible to detect errors in the grouping of substances that have been identified. 

### Preparing data
rm(list = ls())
library(GCalignR)
# aligned data
load("ottensmann-stoffel-hoffman/Supplementary-A/data/beph_aligned.RData")
# linear shifts 
shifts <- beph$Logfile$LinearShift$shift 
beph <- beph$aligned$RT 
# load ms-data and round retention times
beph_ms <- read.table("ottensmann-stoffel-hoffman/Supplementary-A/data/beph_ms.txt",sep = "\t",header = T)
beph_ms[,2:ncol(beph_ms)] <- round(beph_ms[,2:ncol(beph_ms)],digits = 2)
# format as lists
aligned <- as.list(beph[,2:ncol(beph)])
for (i in 1:length(aligned)) {
    aligned[[i]][aligned[[i]] > 0] <- aligned[[i]][aligned[[i]] > 0] - shifts[i]
}
ms <- as.list(beph_ms[,3:ncol(beph_ms)])

## Match peaks to known substances
indices <- list()
# Select the rows that can be linked to known substances for each sample
for (i in 1:length(aligned)) {
 indices <- append(indices,list(which(aligned[[i]] %in% ms[[i]])))   
}
## get all rows
rows <- sort(unique(unlist(indices)))
## update the lists, discard rows that are not informative
fx <- function(x, rows) x[rows]
aligned <- lapply(aligned, fx,rows)
ms <- as.list(beph_ms[,3:ncol(beph_ms)])

merge <- beph[rows,2:ncol(beph)]
rownames(merge) <- 1:length(rownames(merge))

for (i in 1:length(ms)) { # All samples
    for (n in 1:length(ms[[1]])) { # All substances
      if (!is.na(ms[[i]][n])) merge[which(aligned[[i]] == ms[[i]][n]),i] <- as.character(beph_ms$Compounds[[n]])  
    }
}
### Step 2. Define correct row for each substance as the mode, function Mode()

library(MyFunctions)
fy <- function(x,y) { 
    if (any(x %in% y)) { 
        out <- which(x == y)
    } else { 
        out <- NA
    }
    return(out[1])
}   

# set 'correct' rows 
r <- numeric(); x <- 0
for (n in 1:length(beph_ms$Compounds)) {
    r <- append(r, (Mode(apply(merge,2,fy,as.character(beph_ms$Compounds[[n]])))[[1]][[1]]))
    temp <- apply(merge,2,fy,as.character(beph_ms$Compounds[[n]]))
    temp <- temp[!is.na(temp)]
    if (any(temp != r[n])) x <- length(temp[temp != r[n]]) + x
}


error <- x/(nrow(beph_ms) * ncol(beph_ms))
error


