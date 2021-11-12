
#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : FunctionsForGSEA
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : 
#
# createGSEA functions runs the fgsea algorithm of a preranked set against a msigDB collection . Ouputs a table
# with results, a barplot with significant enriched pathways and a sticks plots with the top pathways up and down
#
# tibble_with_lists_to_csv converts the table from createGSEA into a usable csv file
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************



#----------Functions

#function to execute GSEA
createGSEA<-function(statsLab, paths, genetableused ,transition){
  
  library(stats)
  library(fgsea)
  
  library(tidyverse)
  library(dplyr)
  
  fgseaRes <- fgsea(pathways=paths, stats=statsLab,minSize=15,maxSize=500)
  
  #tidy it
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  bardata<-subset(fgseaResTidy,fgseaResTidy$padj<=0.05)
  
  
  #barplot
  barname<-paste(transition,"/",genetableused,"/",genetableused,transition,"barplot.pdf", sep="")
  ggplot(bardata, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=(padj<0.05))) +
    coord_flip() +
    labs(x="Pathways", y="Normalized Enrichment Score",
         title= "NES from GSEA") + 
    theme_minimal()
  ggsave(filename = barname, width = 20, height = 20)
  
  #dev.off()
  
  
  #dim(subset(fgseaResTidy,fgseaResTidy$padj<0.05))
  
  
  #sticks_plot
  
  topPathwaysUp <- fgseaRes[ES > 0, ][head(order(padj), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0, ][head(order(padj), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  stickname<-paste(transition,"/",genetableused,"/",genetableused, transition,"stickstop10.pdf", sep="")
  
  pdf(file = stickname,h=10,w=12)
  plotGseaTable(paths[topPathways], statsLab, fgseaRes,gseaParam = 0.5)
  dev.off()
  
  return(fgseaResTidy)
}


#function to convert the fgeasRes list to dataframe to csv file
tibble_with_lists_to_csv <- function(x, file_path_name) {
  x$leadingEdge<-sapply(x$leadingEdge, paste, collapse=",")
  write.csv(x, file=file_path_name)
}