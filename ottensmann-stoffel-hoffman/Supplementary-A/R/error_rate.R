error_rate <- function(GCalignObject, Reference = "ottensmann-stoffel-hoffman/Supplementary-A/data/bbim_ms.txt") {

### Internal functions
### -----------------------------------------------------------------
    # Mode of a vector
    Mode <- function(x, na.rm = TRUE) {
        if (isTRUE(na.rm)) {
            ux <- unique(x[!is.na(x)])
            x <- x[!is.na(x)]
        } else {
            ux <- unique(x)
        }
        return(ux[which.max(tabulate(match(x, ux)))])
    }
### -----------------------------------------------------------------

    # Obtain linear shifts that have been applied to retention times during alignment
shifts <- GCalignObject[["Logfile"]][["LinearShift"]][["shift"]] 
    # Get retention times
rt <- GCalignObject[["aligned"]][["RT"]]
    # load ms-data and round retention times
ref <- read.table(Reference,sep = "\t",header = T)
ref[,2:ncol(ref)] <- round(ref[,2:ncol(ref)],digits = 2)
    # format as lists
aligned <- as.list(rt[,2:ncol(rt)])
for (i in 1:length(aligned)) {
    aligned[[i]][aligned[[i]] > 0] <- aligned[[i]][aligned[[i]] > 0] - shifts[i]
}
ms <- as.list(ref[,3:ncol(ref)])

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
ms <- as.list(ref[,3:ncol(ref)])

merge <- rt[rows,2:ncol(rt)]
rownames(merge) <- 1:length(rownames(merge))

for (i in 1:length(ms)) { # All samples
    for (n in 1:length(ms[[1]])) { # All substances
        if (!is.na(ms[[i]][n])) merge[which(aligned[[i]] == ms[[i]][n]),i] <- as.character(ref[["Compounds"]][[n]])  
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
for (n in 1:length(ref[["Compounds"]])) {
    r <- append(r, (Mode(apply(merge,2,fy,as.character(ref[["Compounds"]][[n]])))[[1]][[1]]))
    temp <- apply(merge,2,fy,as.character(ref[["Compounds"]][[n]]))
    temp <- as.vector(temp[!is.na(temp)])
    if(any(temp != r[n])) x <- length(temp[temp != r[n]]) + x
}

error <- x/(nrow(ref) * ncol(ref))
return(error)
}