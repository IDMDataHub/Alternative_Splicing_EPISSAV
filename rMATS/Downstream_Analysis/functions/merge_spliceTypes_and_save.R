merge_spliceTypes_and_save <- function(splice_events, new_suffix,colsOI,indir){ 
  
  ##""" Take modified tables with Replicates-columns & Exon-coordinates 
  ##""" and harmonize columns for all splicing-events
  
  listSpliceTabs <- list()
  for (s in splice_events) {
    splice_tab <- read.delim(paste0(indir, s, new_suffix),header=T,as.is=T)
    listSpliceTabs[[s]] <-  splice_tab[,colsOI]
  }
  
  ##---SAVE ----##
  merged_splice_events <- do.call(rbind, listSpliceTabs)
 
  return(merged_splice_events)
  
}