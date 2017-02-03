###########################################
### Functions to simulate chromatograms ###
###########################################

add_peak_error <- function(chroma,p,rt_col_name,conc_col_name,distr = c(-.02,-.01,.01,.02)) {
### add errors to a proportion ("p") of peaks of one sample ("chroma")
    
### 1.) Gaussian distribution of errors
error <- sample(x = distr, 
         replace = T,prob = c(0.15,0.35,0.35,0.15),
         size = round(p * length(chroma[[rt_col_name]][!is.na(chroma[[rt_col_name]])]))) 

### 2.) Select peaks to manipulate
index <- sample(x = 1:length(chroma[[rt_col_name]][!is.na(chroma[[rt_col_name]])]), 
                  size = length(error), 
                  prob = chroma[[conc_col_name]][!is.na(chroma[[conc_col_name]])]/sum(chroma[[conc_col_name]][!is.na(chroma[[conc_col_name]])])) 

### 3.) Get retention times 
retention_time <- chroma[[rt_col_name]][!is.na(chroma[[rt_col_name]])]

### 4.) Define manipulated retention times
retention_time[index] <- retention_time[index] + error 

### 5.) Check for artefacts (i.e. duplicated retention times)
if (any(duplicated(retention_time))) { 
conflict_rt <- which(duplicated(retention_time) == "TRUE") # index of first rt
conflict_rt <- unlist(lapply(conflict_rt, function(x) { # get rt 
   if (!(x %in% index)) { 
    if ((x + 1) %in% index) x <- x + 1
    if ((x - 1) %in% index) x <- x - 1 
  } 
  return(x)
}))
error[conflict_rt] <- error[conflict_rt] * -1 # solve by error *-1
}
 
### 6.) Do the manipulations
  chroma[[rt_col_name]][index] <- chroma[[rt_col_name]][index] + error
  
### 7.) Save errors for subsequent use
x <- numeric(length = length(retention_time))
x[index] <- error
  return(list(chroma = chroma, error = x))
}

correct_rt <- function(optimal_chroma,manipulated_chroma,conc_col_name,rt_col_name) {
# get the original retention times for indexing. Keep positions

    # loop through all samples
for (s in 1:as.numeric(manipulated_chroma[["Logfile"]][["Input"]][["Samples"]])) { 
    area1 <- optimal_chroma[["aligned"]][[conc_col_name]][[s + 1]]  # Optimal aligned
    area1 <- area1[area1 > 0 & !is.na(area1)] # Do not take anything empty
    rt1 <- as.vector(optimal_chroma[["aligned"]][[rt_col_name]][[s + 1]])       
    rt1 <- rt1[rt1 > 0 & !is.na(rt1)]
    
    area2 <- manipulated_chroma[["aligned"]][[conc_col_name]][[s + 1]]    # Manipulated
    area2 <- area2[area2 > 0 & !is.na(area2)]
    rt2 <- as.vector(manipulated_chroma[["aligned"]][[rt_col_name]][[s + 1]])
    rt2 <- rt2[rt2 > 0 & !is.na(rt2)] 
    index <- match(rt2,as.vector(manipulated_chroma[["aligned"]][[rt_col_name]][[s + 1]]),nomatch = 0)
    
    # warning. Need to think what it means!
    if (sum(area2 - area1) != 0) print(paste0("S",s,": Order of Peaks differ between alignments!"))
    area1 <- area1[(match(area2,area1,nomatch = 0))]              
    rt1 <- rt1[(match(area2,area1,nomatch = 0))]
    rt_corrected <- rt1
    
    # change retention times
    manipulated_chroma[["aligned"]][["rt"]][[s + 1]][index] <- rt_corrected
  }
  temp <- rt_to_matrix(manipulated_chroma)  # create martrix of correct retention times
  manipulated_chroma$heatmap_input$aligned_rts[,2:(ncol(temp) + 1)] <- temp # write rts to heatmap
  manipulated_chroma <- correct_heatmap(manipulated_chroma) # change colnames accordingly for heatmap
  return(manipulated_chroma)
}

correct_heatmap <- function(GCalign){
  rt_mat <- GCalign[["heatmap_input"]][["aligned_rts"]] # rt matrix
  rt_mat[,2:ncol(rt_mat)][rt_mat[,2:ncol(rt_mat)] == 0] <- NA # Exclude first row with individual names
  names(GCalign[["heatmap_input"]][["aligned_rts"]])[-1] <- as.character(rt_Mode(rt_mat))[-1]
  return(GCalign) # return the whole object again
}

rt_Mode <- function(x){
  # find mode per column
    name <- rep(NA,ncol(x))
  for (i in 1:ncol(x)) {
    ux <- unique(x[,i])
    ux <- ux[!is.na(ux) & !(ux == 0)]
    ux <- ux[which.max(tabulate(match(x, ux)))]  
    name[i] <- ux
  }
  return(name)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rt_to_matrix <- function(Chromas,rt_col="rt"){
  # writes the retention times from the listed chromatograms to one matrix
  if (class(Chromas) == "GCalign") {
    rt_mat <- matrix(NA,nrow = length(Chromas[["aligned"]][["rt"]]) - 1,ncol = length(Chromas$aligned$rt$S1))
    for (i in 1:nrow(rt_mat)) {
      rt_mat[i,] <- Chromas$aligned$rt[[i + 1]] 
    }
    
    return(rt_mat)
  }  
}

sim_linear_shift <- function(Chromas,shifts,rt_col_name){
  ## Simulates linear shifts in the gas-chromatography by applying 
  ## shifts to all peaks among a chromatogram. Linear shifts are randomly 
  ## sampled from a user-definied range.
  ## Chromas: List of Chromatograms
  ## shifts: Vector defining the imum & maximum value of allowed linear shifts
  
  ## possible linear shifts
  increments <- seq(from = shifts[1],to = shifts[2],by = 0.01) 
  ## Do replicate the random shitfs, a set.seed was used
  set.seed(300)
  ## for each chromatogram except the first acting as a reference.
  shifts <- sample(x = increments,size = length(Chromas) - 1,replace = TRUE)   
  ## The reference is not shifted, otherwise deviations could be larger than interesting
  shifts <- c(0,shifts)
  
  ## Apply the selected shifts
  for (n in 1:length(Chromas)) {
    Chromas[[n]][[rt_col_name]] <-  Chromas[[n]][[rt_col_name]] + shifts[n]
  }
  ## Save a shift matrix for reference
  shift <- matrix(NA,ncol = 2,nrow = length(shifts))
  shift[,2] <- shifts 
  shift[,1] <- names(Chromas)
  colnames(shift) <- c("Sample","Shift")
  return(out = list(Chromas = Chromas,Shifts = shift))
}

na.remove <- function(df) {
    na <- function(x) all(is.na(x))
    index <- which(apply(df,1,na) == TRUE)
    if (length(index) > 0) df <- df[-index,]
    return(df)
}
    
matrix_append <- function(gc_peak_df, gc_peak_list){
    # Add zeros matrices to fit the dimensions of the largest matrix
    MaxLength <- max(sapply(gc_peak_list,function(x) nrow(x)))
    ToAppend <- MaxLength - nrow(gc_peak_df)
    Cols <- ncol(gc_peak_df)
    Zeros <- matrix(0,nrow = ToAppend,ncol = Cols)
    colnames(Zeros) <- names(gc_peak_df)
    gc_peak_df <- rbind(gc_peak_df[,],Zeros)
    return(gc_peak_df)
}

correct_rt <- function(d1 = aligned, d2 = bfla_shifts, rt_col_name = "RT") {
### maybee some checks for reversal of peaks etc.
### Note, linear shifts introduced by align_chroma are not corrected yet.
    rt1 <- d1[["aligned"]][[rt_col_name]][2:length(d1[["aligned"]][[rt_col_name]])] # transformed rts
    rt2 <- d2[[1]] # applied random errors
    rt3 <- data.frame(x = as.numeric(as.character(as.data.frame(d2[[2]])[[2]]))) # linear shifts
    rownames(rt3) <- names(rt2)
    
    for (n in names(rt1)) {
        rt1[[n]][rt1[[n]] > 0] <- rt1[[n]][rt1[[n]] > 0] + rt2[[n]] + rt3[n,]
    }
    d1[["aligned"]][[rt_col_name]][2:length(d1[["aligned"]][[rt_col_name]])] <- rt1
    return(d1)
}

original_rt <- function(org = NULL, aligned = NULL, rt_col_name = NULL){
    # org = input list of alignment a level zero, true rts
    # aligned = GCalign object
    for (i in 1:(length(aligned[["aligned"]][[rt_col_name]]) - 1)) {
        x <- aligned[["aligned"]][[rt_col_name]][[i + 1]] # treated rts
        y <- org[[i]][[rt_col_name]]
        z <- x
        z[z > 0] <- y[y > 0]
        aligned[["aligned"]][[rt_col_name]][i + 1] <- z
    }
    return(aligned)
}
