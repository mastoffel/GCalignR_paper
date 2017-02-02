# bin
# old functions
# 

remove_zero_area_peaks <- function(x){ # remove zeros that were introduced while aligning peaks
    index <- which(x$area==0)# area = 0 --> Peak is absent
    x <- x[-index,] # skip those entries
    return(x)
}

scent2heatmap <- function(scent_matrix){
    # Makes a scent matrix ready for plotting a heatmap
    scent_matrix <- cbind(row.names(scent_matrix),scent_matrix)
    names(scent_matrix)[1] <- "id"
    out <- list(heatmap_input=list(aligned_rts=scent_matrix,initial_rt=scent_matrix))
    return(out)
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

make_dummy_list <- function(Samples=10,Template=template){
    
    # Creates a bunch of chromatograms by sampling peaks from a template chromatogram,
    # acting as a pool of peaks. 
    # Samples: Defines the number of chromatograms to create
    # Peaks: A vector detering the possible number of peaks per chromatogram
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
    # Peaks: A vector detering the possible number of peaks per chromatogram
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
        # Quantifies the correctness of the alignment. For every peaks the row is detered in which the most 
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
                abs(diff(range(rt_mat[,i][!(rt_mat[,i] == 0)],na.rm = T))) 
        }
        return(mean(ranges,na.rm = TRUE))
    }
    
    peaks_per_row <- function(rt_mat) {
        count <- rep(NA,ncol(rt_mat))
        for(i in 1:ncol(rt_mat)) {
            count[i] <- 
                length(rt_mat[,i][rt_mat[,i]])
        }
        return(count)
    }

#     set.seed(300)
#     indices <- lapply(N, index,Template) # Indices of Peaks for every sample
#     
#     counter <- 1
#     for (i in seq(from = 1,by = 2,to = Samples*2)) {
#         peak_mat[1:length(indices[[counter]]),i] <- Template[[rt_col_name]][indices[[counter]]]
#         peak_mat[1:length(indices[[counter]]),i + 1] <- Template[[conc_col_name]][indices[[counter]]]
#         counter <- counter + 1 
#     }
#     
#     peak_mat <- as.data.frame(peak_mat)
#     indices <- seq(from = 1,by = 2,to = Samples*2)
#     peak_list <- lapply(indices, function(x) peak_mat[,x:(x + 1)])
#     names(peak_list) <- paste0("S",as.character(1:Samples))
#     
#     rename_cols = function(data, var_names){
#         names(data) <- var_names
#         data
#     }
#     
#     peak_list <- lapply(peak_list, rename_cols, c(rt_col_name,conc_col_name))
#     return(peak_list)
# } 

