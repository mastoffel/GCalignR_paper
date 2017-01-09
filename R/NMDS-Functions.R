# MDS Functions, Data wrangling, plotting and statistics


myMetaMDS <- function(scent,factors){
# Does a NMDS approach and combines scent abundance data with factors in a data.frame
scent <- scent[match(row.names(factors),row.names(scent)),] # sort by factors
scent_nmds <- vegan::metaMDS(scent, distance = "bray", k=2, trymax=1000,autotransform=FALSE, expand=FALSE, plot=F)
scent_nmds <- as.data.frame(scent_nmds$points)
scent_nmds <- cbind(scent_nmds,factors)
scent_nmds$family <- as.factor(scent_nmds$family)
return(scent_nmds)  
}

scent_extract <- function(GCalignR_Object,factors,conc_col_name="area",rt_col_name="rt"){
  scent <- log(norm_peaks(GCalignR_Object,conc_col_name,rt_col_name,out = "data.frame")+1) # normalise and log-transform
  scent <- scent[match(row.names(factors),row.names(scent)),] # sort by factors
return(scent)
}

plot_nmds_colony <- function(nmds_df){
  p <- 
    ggplot2::ggplot(data = nmds_df,aes(MDS1,MDS2,color=as.factor(colony))) +
    stat_ellipse(size=2) +
    geom_point(size=7) +
    # geom_text(aes(label=family),size=4,colour="gray25") +
    scale_x_continuous(limits = c(-1.25,1.25)) +
    scale_y_continuous(limits = c(-1.25,1.25)) +
    labs(title ="",
         x = "MDS1", 
         y = "MDS2")+
    theme_bw()+theme(
      plot.title=element_text(face = "bold"),
      axis.title.x = element_text(size = 16,vjust=0.2),
      axis.text.x  = element_blank(), 
      axis.title.y = element_text(size = 16,vjust=0.2),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank())
      # geom_rect(colour="black",fill="khaki1",aes(xmin=1,xmax=1.2,ymin=0.95,ymax=1.25),alpha=0.2) +
      # annotate("text",x=1.1,y=1.1,label=paste(c("Adonis",R2,Pval),collapse = "\n"),col="Black",size=5)
  return(p)
}

plot_nmds_family <- function(nmds_df){
  p <- 
    ggplot2::ggplot(data = nmds_df,aes(MDS1,MDS2,color=as.factor(colony))) +
    # stat_ellipse(size=2) +
    geom_point(size=7) +   
    geom_line(aes(group = family), size = 1,linetype="dotted",color="black") +
    #geom_text(aes(label=family),size=4,colour="black") +
    scale_x_continuous(limits = c(-1.25,1.25)) +
    scale_y_continuous(limits = c(-1.25,1.25)) +
    labs(title ="",
         x = "MDS1", 
         y = "MDS2")+
    theme_bw()+theme(
      plot.title=element_text(face = "bold"),
      axis.title.x = element_text(size = 16,vjust=0.2),
      axis.text.x  = element_blank(), 
      axis.title.y = element_text(size = 16,vjust=0.2),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="none",
      panel.grid = element_blank())
  # geom_rect(colour="black",fill="khaki1",aes(xmin=1,xmax=1.2,ymin=0.95,ymax=1.25),alpha=0.2) +
  # annotate("text",x=1.1,y=1.1,label=paste(c("Adonis",R2,Pval),collapse = "\n"),col="Black",size=5)
  return(p)
}


# myAdonis <- function(df){
#   scent <- df  
#   stats <- adonis(scent~factors$family,permutations=999,method="bray")
#   
# }

adonis_family <- function(df){
  scent <- df  
  stats <- adonis(scent~factors$family,permutations=999,method="bray",strata = factors$colony)
  
}

adonis_colony <- function(df){
  scent <- df  
  stats <- adonis(scent~factors$colony,permutations=999,method="bray")
  
}