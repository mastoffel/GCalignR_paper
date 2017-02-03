error_rate <- function(GCalignObject, Reference = NULL, rt_col_name = "RT", linshift = TRUE) {

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

    
    # Get retention times
rt <- GCalignObject[["aligned"]][[rt_col_name]]
    # load ms-data and round retention times
ref <- read.table(Reference,sep = "\t",header = T)
ref[,2:ncol(ref)] <- round(ref[,2:ncol(ref)],digits = 2)
    # format as lists
aligned <- as.list(round(rt[,2:ncol(rt)],2))
if (isTRUE(linshift)) {
    # Obtain linear shifts that have been applied to retention times during alignment
    shifts <- GCalignObject[["Logfile"]][["LinearShift"]][["shift"]] 
for (i in 1:length(aligned)) {
    aligned[[i]][aligned[[i]] > 0] <- aligned[[i]][aligned[[i]] > 0] - shifts[i]
}
    aligned <- lapply(aligned,round,digits = 2)
} 
ms <- as.list(ref[,3:ncol(ref)])

## Match peaks to known substances
indices <- list()
# Select the rows that can be linked to known substances for each sample
for (i in 1:length(aligned)) {
    indices <- append(indices,list(which(aligned[[i]][aligned[[i]] > 0] %in% ms[[i]][ms[[i]] > 0])))  
    t1 <- length(which(aligned[[i]][aligned[[i]] > 0] %in% ms[[i]][ms[[i]] > 0]))
    t2 <- length(ms[[i]][!is.na(ms[[i]]) & ms[[i]] > 0])
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
    if (any(temp != r[n])) x <- length(temp[temp != r[n]]) + x
}

temp <- unlist((ref[1:nrow(ref),3:ncol(ref)]))
N <- length(temp[!is.na(temp) & temp > 0])
N1 <- length(unlist(indices))
if (N > N1) warning("Not all retention times were matched")
# check
error <- x/N
return(error)
}
