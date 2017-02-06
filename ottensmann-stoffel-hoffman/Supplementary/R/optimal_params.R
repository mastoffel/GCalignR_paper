## Determine optimal alignment paramters by minimising errors
## returns aligned data and a data frame listing the paramter combination
optimal_params <- function(data,rt_col_name,conc_col_name,max_diff_peak2mean,min_diff_peak2peak) {
    library(GCalignR)
    out <- list()
    names <- character()
    params <- matrix(NA, nrow = length(max_diff_peak2mean)*length(min_diff_peak2peak),ncol = 2)
    counter <- 1
    for (i in max_diff_peak2mean) {
        for (n in min_diff_peak2peak) {
            out <- append(out,list(align_chromatograms(data = data,rt_col_name = rt_col_name,max_diff_peak2mean = i,min_diff_peak2peak = n,max_linear_shift = 0.01)))
            names <- c(names,paste0("d_","p2m_",as.character(i),"_p2p_",as.character(n)))
            params[counter,1] <- i
            params[counter,2] <- n
            counter <- counter + 1
                }
    }
    names(out) <- names
    params <- data.frame(p2m = params[,1],p2p = params[,2])
    return(list(out,params))
}

    
