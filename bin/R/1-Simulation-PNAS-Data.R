## Do simulations based on the PNAS Data
## Assumptions:
## The Alignment performed with GCalignR on the raw resembles the "true" pattern 
## i.e. a strong location effect of the colony, and no dispersion effect of data!
## 
## DO NOT run this script, as it takes long time. Additionally it might lack some 
## functionallity, because of improved GCalignR code

# A: Transform the aligned Data to the input format of GCalignR
#    Output is a list called Chromas with data.frame for Inds
# --------------------------------------------------------------

rm(list = ls())
## Load packages & functions
library(GCalignR); source("R/ChromaSimFunctions.R")
## Read GCalignR output
load("data/PNAS/afs_peaks_aligned_100816.RData") 
## Inspect the alignment accuracy with a heatmap
gc_heatmap(object = afs,threshold = 0.02)

## Extract the area of peaks for all chromatograms
area <- as.data.frame(as.list(afs$aligned$area[2:length(afs$aligned$area)]))
## Calculate the mean retention time for every substance, these are used for every sample!
mean_rt <- round(afs$aligned$rt$mean_RT,digits = 2)

## build a list of samples
Chromas <- lapply(1:ncol(area), function(x) as.data.frame(cbind(mean_rt,area[,x]))) 
## rename every sample accordingly
names(Chromas) <- names(area) 
rename_cols <- function(data, var_names){
  names(data) <- var_names
  return(data) }

## name the variables "rt" and "area"
Chromas <- lapply(Chromas, rename_cols, c("rt","area"))
## remove zero area peaks
Chromas <- lapply(Chromas, remove_zero_area_peaks) 

### B. Apply Linear- Shifts, these are applied randomly, but are saved to trace back
## The range of shifts is restricted between -0.03 and 0.03
out <- sim_linear_shift(Chromas = Chromas,shifts = c(-0.03,0.03))
## Get List of systematically shifted retention times
ChromasLinShift <- out$Chromas
Applied_Shifts <- out$Shifts

## clean workspace
rm(list = setdiff(ls(), c("ChromasLinShift","Chromas","Applied_Shifts")))
## Load stuff again
source("R/ChromaSimFunctions.R"); library(GCalignR)

### Start the simulations
### For speed reasons, the simulations was divided in ten parts, that were run in parallel

part <- 1 # 1:10 were all done!

## Define proportions of random errors: The value determines how many peaks (closest to the desired proportion) are treated by adding random errors in the range -0.02:0.02, without zero
peak_quality <- rep(seq(from = 0,to = 1,by = 0.1),each = 10) 

### C: Align the list of untreated chromatograms & the manipulated list 

## OptAlign is a reference, no systematic shifts, no random shifts
OptAlign <- align_chromatograms(Chromas,
                                rt_col_name = "rt",
                                conc_col_name = "area",
                                reference = "P39") 

## for every level of peak quality (0 to 1) and its replicates, 1110 in total
for (i in 1:length(peak_quality)) { 
  
  ## Add random errors
  ChromasGaussErr <- lapply(ChromasLinShift,add_peak_error,p = peak_quality[i]) 

  ## Align the manipulated chromatogram data
  SimAlign <- align_chromatograms(ChromasGaussErr,
                                  rt_col_name = "rt",
                                  conc_col_name = "area",
                                  reference = "P39")   
## Write the GCalingR output to a file
assign(paste0("sim_",as.character(i)),SimAlign)
}


SimAlign <- list(sim_1 = sim_1,
                 sim_2 = sim_2,
                 sim_3 = sim_3,
                 sim_4 = sim_4,
                 sim_5 = sim_5,
                 sim_6 = sim_6,
                 sim_7 = sim_7,
                 sim_8 = sim_8,
                 sim_9 = sim_9,
                 sim_10 = sim_10,
                 sim_11 = sim_11,
                 sim_12 = sim_12,
                 sim_13 = sim_13,
                 sim_14 = sim_14,
                 sim_15 = sim_15,
                 sim_16 = sim_16,
                 sim_17 = sim_17,
                 sim_18 = sim_18,
                 sim_19 = sim_19,
                 sim_20 = sim_20,
                 sim_21 = sim_21,
                 sim_22 = sim_22,
                 sim_23 = sim_23,
                 sim_24 = sim_24,
                 sim_25 = sim_25,
                 sim_26 = sim_26,
                 sim_27 = sim_27,
                 sim_28 = sim_28,
                 sim_29 = sim_29,
                 sim_30 = sim_30,
                 sim_31 = sim_31,
                 sim_32 = sim_32,
                 sim_33 = sim_33,
                 sim_34 = sim_34,
                 sim_35 = sim_35,
                 sim_36 = sim_36,
                 sim_37 = sim_37,
                 sim_38 = sim_38,
                 sim_39 = sim_39,
                 sim_40 = sim_40,
                 sim_41 = sim_41,
                 sim_42 = sim_42,
                 sim_43 = sim_43,
                 sim_44 = sim_44,
                 sim_45 = sim_45,
                 sim_46 = sim_46,
                 sim_47 = sim_47,
                 sim_48 = sim_48,
                 sim_49 = sim_49,
                 sim_50 = sim_50,
                 sim_51 = sim_51,
                 sim_52 = sim_52,
                 sim_53 = sim_53,
                 sim_54 = sim_54,
                 sim_55 = sim_55,
                 sim_56 = sim_56,
                 sim_57 = sim_57,
                 sim_58 = sim_58,
                 sim_59 = sim_59,
                 sim_60 = sim_60,
                 sim_61 = sim_61,
                 sim_62 = sim_62,
                 sim_63 = sim_63,
                 sim_64 = sim_64,
                 sim_65 = sim_65,
                 sim_66 = sim_66,
                 sim_67 = sim_67,
                 sim_68 = sim_68,
                 sim_69 = sim_69,
                 sim_70 = sim_70,
                 sim_71 = sim_71,
                 sim_72 = sim_72,
                 sim_73 = sim_73,
                 sim_74 = sim_74,
                 sim_75 = sim_75,
                 sim_76 = sim_76,
                 sim_77 = sim_77,
                 sim_78 = sim_78,
                 sim_79 = sim_79,
                 sim_80 = sim_80,
                 sim_81 = sim_81,
                 sim_82 = sim_82,
                 sim_83 = sim_83,
                 sim_84 = sim_84,
                 sim_85 = sim_85,
                 sim_86 = sim_86,
                 sim_87 = sim_87,
                 sim_88 = sim_88,
                 sim_89 = sim_89,
                 sim_90 = sim_90,
                 sim_91 = sim_91,
                 sim_92 = sim_92,
                 sim_93 = sim_93,
                 sim_94 = sim_94,
                 sim_95 = sim_95,
                 sim_96 = sim_96,
                 sim_97 = sim_97,
                 sim_98 = sim_98,
                 sim_99 = sim_99,
                 sim_100 = sim_100,
                 sim_101 = sim_101,
                 sim_102 = sim_102,
                 sim_103 = sim_103,
                 sim_104 = sim_104,
                 sim_105 = sim_105,
                 sim_106 = sim_106,
                 sim_107 = sim_107,
                 sim_108 = sim_108,
                 sim_109 = sim_109,
                 sim_110 = sim_110) 

### E:) Save Output           

SimData < -list(OptAlign = OptAlign,
              SimAlign = SimAlign,Linshift = Applied_Shifts,peak_quality = peak_quality) 
save(x = SimData,file = paste0("data/",Sys.Date(),"-",as.character(part),"-PNAS",".RData"))   
 
