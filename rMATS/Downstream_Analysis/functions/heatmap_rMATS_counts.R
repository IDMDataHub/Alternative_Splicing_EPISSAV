#!/usr/bin/env R 

# source("manip_rmats_output_star.R")

plot_heatmap_signif_splice <- function(ase, fdr.co, psi.co, pair, plot_name) {
  
    ase.filt <- subset(ase, (FDR < fdr.co &  abs(IncLevelDifference) > psi.co) )
    
    
    ###  z.gene =  (E(g) - Mean(g.replicates) ) / M.A.D (g.replicates), mad=median.absolute.deviation
    x <- ase.filt[,8:13]
    s1 = strsplit(pair, "_vs_")[[1]][1]
    s2 = strsplit(pair, "_vs_")[[1]][2]
    colnames(x) <- c( paste(s1,"_", c("r1","r2","r3"),sep=""),
                      paste(s2,"_",c("r1","r2","r3"),sep=""))
    
    
    cal_z_score <- function(x){
      (x - mean(x)) / sd(x)
    }
    x_norm <- t(apply(x, 1, cal_z_score))
    # x_norm <- x_norm[-which(apply(x_norm, 2, is.na)) ,]
    
    #Rowv=as.dendrogram(hcl_row),tree=ct))
    library(gplots)
    library(RColorBrewer)
    
   
    heatmap.2(x_norm, 
              Rowv=TRUE,
              Colv=F,
              dendrogram ='row',
              scale="row", na.rm=TRUE,
              trace='none',
              labCol=colnames(x),
              col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
              margins=c(9,2),
              cexCol=1.3,
              srtCol=45,
              labRow="",
              offsetCol=0,
              key=TRUE,
              keysize = 1.5,
              density.info="density",
              key.title="Color range",
              key.xlab="norm.counts",
              key.ylab="Counts",
              key.par=list(mar=c(3,4,3,1)),
              notecex=3, 
              main=paste(nrow(x_norm),"signif. Diff/ial splicing events\nin ",pair))
    #print(hmp)
   
}


## FILTER NON-SIGNIFICANT 

fdr.co <- 0.05
#psi.co <- 0.15
psi.co <- 0.20

pair="Saline_vs_MBNLdecoy"
plot_name="Plots/Heatmap_SignifAltSplicing_Encode_Saline_vs_MBNLdecoy.png"

# png(plot_name)
plot_heatmap_signif_splice(ase.saline.filt, fdr.co, psi.co, pair, plot_name)
#dev.off()

## 18 events when psi.co = 0.15
## 14 events when psi.co = 0.20


####----- PIE-CHART of SIGNIF-SPLICED EVENTS and Others -------####
#### Show that differential splicing events in MBNLdecoy compaired to Saline are very few 

create_pieChart_Splicing <- function(ase,pair,signif.co,psi.co) {
  ##""" Test with : ase = ase.saline; pair="Saline_vs_MBNLdecoy" ; signif.co =0.05 ; psi.co=0.20 """
  
  library(ggplot2)

  ## Count events Signif. splicing : 
  ## count InclDifference(dPSI) = mean(Saline) -  mean(MBNLdecoy)
  ## INCLUDED : PSI < -0.15 
  ## EXCLUDED : PSI > 0.15
  
  signif <- nrow(subset(ase, (FDR < signif.co  &  abs(IncLevelDifference) >= psi.co )))
  all_splice <-  nrow(ase)
  
  splice_counts <- data.frame(category=c("All", "Signif & highPSI"),
                              count=c(all_splice,signif))
  
  ## calculate percentage of each count : 
  splice_counts$fraction = round(splice_counts$count / sum(splice_counts$count),2) * 100
  
  # Compute the cumulative percentages (top of each rectangle)
  splice_counts$ymax = cumsum(splice_counts$fraction)
  
  # Compute the bottom of each rectangle
  splice_counts$ymin = c(0, head(splice_counts$ymax, n=-1))
  
  # Compute label position
  # pos <- (splice_counts$ymax + splice_counts$ymin) / 2 # = 49.0 98.5 99.5
  # 98 & 99 are very close, I ll replace them with 97 and 10, left and right of the end of circle.
  splice_counts$labelPosition <- c(50.0, 97.0)
  
  # Compute a good label
  splice_counts$label <- paste0(splice_counts$category,"\n", splice_counts$fraction, " %")
  
  # Make the plot
  dir.create("Plots/",showWarnings = F)
  #png(paste0("Plots/Donut_Plot_Signif_SplicingEvents_",pair,".png"))
    ggplot(splice_counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect() +
      geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") + 
      xlim(c(2, 4)) + 
      theme_void()
 #dev.off()
    
}

## call : 
create_pieChart_Splicing(ase.saline.star, pair="Saline_vs_MBNLdecoy",signif.co=0.05, psi.co=0.20)




