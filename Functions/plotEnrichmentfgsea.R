

#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : plotEnrichmentfgsea

#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : 
# 
#  plotEnrichmentfgsea improves the original enrichment plot from the fgsea package
#
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************



#modify plot


plotEnrichmentfgsea <- function(pathway, stats,
                                gseaParam=1,
                                ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="green", size=1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed", size =1) +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed", size =1) +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="green") + theme_bw(base_size = 24) +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(x="Rank", y="Enrichment Score")
  g
}