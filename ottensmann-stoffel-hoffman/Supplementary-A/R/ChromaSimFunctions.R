###########################################
### Functions to simulate chromatograms ###
###########################################

### add errors to a proportion ("p") of peaks of one sample ("chroma")
add_peak_error <- function(chroma,p,rt_col_name,conc_col_name) {

### 1.) Gaussian distribution of errors
error <- sample(x = c(-0.02,-0.01,0.01,0.02), 
         replace = T,prob = c(0.15,0.35,0.35,0.15),
         size = round(p * length(chroma[[rt_col_name]][!is.na(chroma[[rt_col_name]])]))) 

### 2.) Select peaks to manipulate
index <- sample(x = 1:length(chroma[[rt_col_name]][!is.na(chroma[[rt_col_name]])]), 
                  size = length(error), 
                  prob = chroma[[conc_col_name]][!is.na(chroma[[conc_col_name]])]/sum(chroma[[conc_col_name]][!is.na(chroma[[conc_col_name]])])) 

### 3.) Get retention times 
retention_time <- chroma[[rt_col_name]] 

### 4.) Define manipulated retention times
retention_time[index] <- chroma[[rt_col_name]][index] + error 

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
  return(chroma)
}

clean_chroma <- function(ChromaList,ChromaRef,conc_col_name, rt_col_name) { 
  # remove duplicates in chroma, by the help of area indices in ChromaGauss   
  for (i in 1:length(ChromaList)) {
    areas <- setdiff(ChromaList[[i]][[conc_col_name]], ChromaRef[[i]][[conc_col_name]])
    indices <- which(ChromaList[[i]][[conc_col_name]] %in% areas)
    ChromaList[[i]][[conc_col_name]][indices] <- NA
    ChromaList[[i]][[rt_col_name]][indices] <- NA
  }
  ChromaList <- lapply(ChromaList,function(x) x[with(x, order(conc_col_name, rt_col_name)),]) # reorder
}

#######################################################################  
# Get true rt & area and their positions in the aligned chromatograms #
# without any manipulation (No linear shift or gaussian error added)  #
#######################################################################

align_correctness <- function(aligned_chroma,optimal_chroma,conc_col_name,rt_col_name){
  quality <- matrix(NA,nrow = as.numeric(aligned_chroma[["Logfile"]][["Input"]][["Samples"]]), ncol = 1)
  for (s in 1:nrow(quality)) { # for every chromatoram
    
    x1 <- as.vector(optimal_chroma[["aligned"]][[rt_col_name]][[s + 1]]) # optimally aligned
    peaks1 <- x1[x1 > 0] # extract all true peaks 
    y1 <- optimal_chroma[["aligned"]][[conc_col_name]][[s + 1]]
    index1 <- match(peaks1,x1)
    area1 <- y1[index1] 
    
    x2 <- as.vector(aligned_chroma[["aligned"]][[rt_col_name]][[s + 1]])     # Manipulated
    peaks2 <- x2[x2 > 0] # extract all true peaks 
    y2 <- aligned_chroma[["aligned"]][["area"]][[s + 1]]
    index2 <- match(peaks2,x2)
    area2 <- y2[index2] 
    
    correct <- index2 - index1
    quality[s] <- length(correct[correct == 0])/length(peaks1) 
  }
  return(quality)
} 

correct_rt <- function(optimal_chroma,manipulated_chroma,conc_col_name,rt_col_name) {
  # Dertermine the cumulative error that was added to every peak by linear- and non linear
  # manipulations of the retention time and correct this in the final alignment. Peaks remain
  # where they have been sorted in!
  # optimal_Chroma: Alignment of the optimal chromatograms. Retention times are unchanged
  # manipulated_chroma: Alignment of manipulated chromatograms. Retention times have been changed randomly
  
  for (s in 1:as.numeric(manipulated_chroma[["Logfile"]][["Input"]][["Samples"]])) { 
    
    area1 <- optimal_chroma[["aligned"]][[conc_col_name]][[s + 1]]       # Optimal aligned
    area1 <- area1[area1 > 0 & !is.na(area1)] # Do not take anything empty
    peaks1 <- as.vector(optimal_chroma[["aligned"]][[rt_col_name]][[s + 1]])       
    peaks1 <- peaks1[peaks1 > 0 & !is.na(peaks1)]
    
    area2 <- manipulated_chroma[["aligned"]][[conc_col_name]][[s + 1]]    # Manipulated
    area2 <- area2[area2 > 0 & !is.na(area2)]
    peaks2 <- as.vector(manipulated_chroma[["aligned"]][[rt_col_name]][[s + 1]])
    peaks2 <- peaks2[peaks2 > 0 & !is.na(peaks2)] 
    index <- match(peaks2,as.vector(manipulated_chroma[["aligned"]][[rt_col_name]][[s + 1]]),nomatch = 0)
    
    if (sum(area2-area1)!=0) print(paste0("S",s,": Order of Peaks differ between alignments!"))
    
    
    #######################
    # Backtransform peaks #
    #######################
    area1 <- area1[(match(area2,area1,nomatch = 0))]              # Sort with respect to area2
    peaks1 <- peaks1[(match(area2,area1,nomatch = 0))]
    
    peaks_corrected <- peaks1
    
    #####################################################
    # Overwrite retention times with the correct values #
    #####################################################
    manipulated_chroma[["aligned"]][["rt"]][[s + 1]][index] <- peaks_corrected
  }
  temp <- rt_to_matrix(manipulated_chroma)  # create martrix of correct retention times
  manipulated_chroma$heatmap_input$aligned_rts[,2:(ncol(temp)+1)] <- temp # write rts to heatmap
  manipulated_chroma <- correct_heatmap(manipulated_chroma) # change colnames accordingly for heatmap
  return(manipulated_chroma)
}

correct_heatmap <- function(GCalign){
  # correct the colnames of the heatmap to allow the calculation of classification scores
  # For determining the correct colname with rt_Mode all values 0 are coded as NA
  
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

GCalign_to_matrix <- function(GCalign,rt_col_name){
  L <- length(GCalign[["aligned"]][[rt_col_name]]) - 1
  rt_mat <- matrix(data = NA,nrow = L,ncol = length(GCalign[["aligned"]][[rt_col_name]][[1]]))
  for (i in 1:L) {
    rt_mat[i,] <- GCalign[["aligned"]][[rt_col_name]][[i + 1]]
  }
  return(rt_mat)
}

make_chroma_list <- function(Samples,Peaks,Template = template,rt_col_name,conc_col_name){
  # Creates a bunch of chromatograms by sampling peaks from a template chromatogram,
  # acting as a pool of peaks. 
  # Samples: Defines the number of chromatograms to create
  # Peaks: A vector determining the possible number of peaks per chromatogram
  # Template: A list containing retention times "rt" and "area"
  # Set.seed ensures to replicate the pseudorandom chromatograms again, when picking the 
  # very same settings
  
  peak_mat <- matrix(data = NA,nrow = max(Peaks),ncol = Samples*2)
  # In peak_mat two columns represent rt & area of a sample. Samples are concatenated
  # column-wise
  
  set.seed(300)
  N <- sample(Peaks,size = Samples,replace = T) # Sizes of Chromatograms
  
  index <- function(N,Template){ # Function for picking peaks
  indices <- sort(sample(x = 1:length(Template[[rt_col_name]]),size = N)) 
  }
  
  set.seed(300)
  indices <- lapply(N, index,Template) # Indices of Peaks for every sample
  
  counter <- 1
  for (i in seq(from = 1,by = 2,to = Samples*2)) {
    peak_mat[1:length(indices[[counter]]),i] <- Template[[rt_col_name]][indices[[counter]]]
    peak_mat[1:length(indices[[counter]]),i + 1] <- Template[[conc_col_name]][indices[[counter]]]
    counter <- counter + 1 
  }
  
  peak_mat <- as.data.frame(peak_mat)
  indices <- seq(from = 1,by = 2,to = Samples*2)
  peak_list <- lapply(indices, function(x) peak_mat[,x:(x + 1)])
  names(peak_list) <- paste0("S",as.character(1:Samples))
  
  rename_cols = function(data, var_names){
    names(data) <- var_names
    data
  }
  
  peak_list <- lapply(peak_list, rename_cols, c(rt_col_name,conc_col_name))
  return(peak_list)
} 

peak_counter <- function(Chromas){
  # Peakcounter estimates how many peaks are called within all chromatograms
  
  rt <- numeric(0)
  for (i in 1:length(Chromas)) {
    rt <- c(rt,Chromas[[i]][["rt"]])
  }
  rt <- unique(rt[!(is.na(rt)) & rt != 0])
  return(length(rt))
}

peak_lister <- function(Chromas,rt_col_name){
  rt <- numeric(0)
  for (i in 1:length(Chromas)) {
    rt <- c(rt,Chromas[[i]][[rt_col_name]])
  }
  rt <- unique(rt[!(is.na(rt)) & rt != 0])
  return(sort(rt))
}

quantify_align <- function(GCalign){
  # Set of functions for quantifying how well the alignment sorted peaks.
  # Therefore for every peak the best fitting row is defined as the row containing
  # the greates proportion of this peak. The quality score is then calculated as the ratio
  # of peaks in that row to the number of samples carrying the peak
  # E.g. 4/6 --> 75 % correctly aligned
  # Quantifies the correctness of the alignment. For every peaks the row is determined in which the most 
  # carrying samples have the respective peak. Therefore a score indicating the proportion of correctly aligned
  # peaks is estimated
  
#   rt <- unique(names(GCalign$heatmap_input$aligned_rts)[-1])
  rt_mat <- GCalign$heatmap_input$aligned_rts[,-1]
  rt <- unique(unlist(rt_mat[!(rt_mat==0)&!is.na(rt_mat)]))
  score <- unlist(lapply(rt, peak_search,rt_mat))
  return(score)
}

peak_search <- function(rt,rt_mat){ # lapply in the end!  
  hits <- as.data.frame(which(rt_mat==rt,arr.ind = T)) # where are the peaks located ?
  N <- nrow(hits)
  TrueCol <- Mode(hits$col) # TrueCol is defined as the column contaning most of the peaks
  score <- length(hits$col[hits$col==TrueCol])/N # percentage grouped correctly
  return(score)  
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rel_var <- function(GCalign){
  # calculates the variation coefficient
  rt_mat <- GCalign_to_matrix(GCalign = GCalign) # get the retention times
  varcoeff <- (rep(NA,ncol(rt_mat)))
  for (i in 1:length(varcoeff)){
    m <- mean(rt_mat[,i][!(rt_mat[,i]==0)],na.rm = TRUE)
    s <- sd(rt_mat[,i][!(rt_mat[,i]==0)],na.rm = TRUE)
    varcoeff[i] <- round((s/m),digits = 2) # coefficient of variation
  }
  return(round(max(varcoeff,na.rm = TRUE),digits = 2))
}

rt_range <- function(GCalign){
  # Calculates the mean range of retention times per row.
  # Singles Peaks, only present in one row are skipped to 
  # avoid inflation of zeros
  rt_mat <- GCalign_to_matrix(GCalign = GCalign) # get the retention times
  ranges <- (rep(NA,ncol(rt_mat)))
  count <- peaks_per_row(rt_mat)
  for (i in 1:length(ranges)){
    ranges[i] <-
      abs(diff(range(rt_mat[,i][!(rt_mat[,i]==0)],na.rm = T))) 
  }
  return(mean(ranges,na.rm = TRUE))
}

peaks_per_row <- function(rt_mat){
  count <- rep(NA,ncol(rt_mat))
  for(i in 1:ncol(rt_mat)){
    count[i] <- 
      length(rt_mat[,i][rt_mat[,i]])
  }
  return(count)
}

rt_to_matrix <- function(Chromas,rt_col="rt"){
  # writes the retention times from the listed chromatograms to one matrix
  if(class(Chromas)=="GCalign"){
    rt_mat <- matrix(NA,nrow =length(Chromas[["aligned"]][["rt"]])-1,ncol = length(Chromas$aligned$rt$S1))
    for(i in 1:nrow(rt_mat)){
      rt_mat[i,] <- Chromas$aligned$rt[[i+1]] 
    }
    
    return(rt_mat)
  }  
}

sim_linear_shift <- function(Chromas,shifts){
  ## Simulates linear shifts in the gas-chromatography by applying 
  ## shifts to all peaks among a chromatogram. Linear shifts are randomly 
  ## sampled from a user-definied range.
  ## Chromas: List of Chromatograms
  ## shifts: Vector defining the minimum & maximum value of allowed linear shifts
  
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
    Chromas[[n]][["rt"]] <-  Chromas[[n]][["rt"]] + shifts[n]
  }
  ## Save a shift matrix for reference
  shift <- matrix(NA,ncol = 2,nrow = length(shifts))
  shift[,2] <- shifts 
  shift[,1] <- names(Chromas)
  colnames(shift) <- c("Sample","Shift")
  return(out = list(Chromas = Chromas,Shifts = shift))
}

make_dummy_list <- function(Samples=10,Template=template){
  
  # Creates a bunch of chromatograms by sampling peaks from a template chromatogram,
  # acting as a pool of peaks. 
  # Samples: Defines the number of chromatograms to create
  # Peaks: A vector determining the possible number of peaks per chromatogram
  # Template: A list containing retention times "rt" and areas "area"
  
  peak_mat <- matrix(data = NA,nrow = 10,ncol = Samples*2)
  # In peak_mat two columns represent rt & area of a sample. Samples are concatenated
  # column-wise
  
  counter <- 1
  for(i in seq(from=1,by=2,to=Samples*2)){
    indices <- counter:(counter+9)  
    peak_mat[1:length(indices),i] <- Template[["rt"]][indices]
    peak_mat[1:length(indices),i+1] <- Template[["area"]][indices]
    counter <- counter+10
  }
  
  peak_mat <- as.data.frame(peak_mat)
  indices <- seq(from=1,by=2,to=Samples*2)
  peak_list <- lapply(indices, function(x) peak_mat[,x:(x+1)])
  names(peak_list) <- paste0("S",as.character(1:Samples))
  
  rename_cols = function(data, var_names){
    names(data) <- var_names
    data
  }
  
  peak_list <- lapply(peak_list, rename_cols, c("rt","area"))
  return(peak_list)
} 

MinMax <- function(x,y){ # Estimate the range of retention times per substance, they should be no overlapp
  # x: Matrix of retention times
  # y: Means per row (i.e. name of the substance)
  temp <- matrix(NA,1,4)
  colnames(temp) <- c("min","max","rt","row")
  data <- temp[0,]
  for(i in 1:ncol(x)){
    data<-rbind(data,cbind(min(x[,i][x[,i]>0],na.rm = TRUE),max(x[,i][x[,i]>0],na.rm = TRUE),y[i],i))
  }

  return(as.data.frame(data))
}

PeakOverlap <- function(rt_min_max,threshold=0.02){
  rt_min_max["width"] <- rt_min_max["max"] - rt_min_max["min"] # correct
  rt_min_max["distance"] <- 0 # preallocate
  for(n in 2:nrow(rt_min_max)){
    rt_min_max[n,"distance"] <- rt_min_max[n,"min"] - rt_min_max[n-1,"max"]
  }
  rt_min_max["distance"][ rt_min_max["distance"]>threshold] <- threshold
  rt_min_max["offset"] <- 0 
  for(n in 2:nrow(rt_min_max)){
  rt_min_max[n,"offset"] <- rt_min_max[n,"distance"] + rt_min_max[n-1,"width"] 
  }
  rt_min_max["cumul"] <- cumsum(rt_min_max["distance"]+rt_min_max["offset"])  

  for(n in 2:nrow(rt_min_max)){
    rt_min_max[n,"offset"] <- rt_min_max[n,"distance"] + rt_min_max[n-1,"width"] + rt_min_max[n-1,"cumul"]
  }
  return(rt_min_max)
 }

scent2heatmap <- function(scent_matrix){
  # Makes a scent matrix ready for plotting a heatmap
  scent_matrix <- cbind(row.names(scent_matrix),scent_matrix)
  names(scent_matrix)[1] <- "id"
  out <- list(heatmap_input=list(aligned_rts=scent_matrix,initial_rt=scent_matrix))
  return(out)
}


remove_zero_area_peaks <- function(x){ # remove zeros that were introduced while aligning peaks
  index <- which(x$area==0)# area = 0 --> Peak is absent
  x <- x[-index,] # skip those entries
  return(x)
}
