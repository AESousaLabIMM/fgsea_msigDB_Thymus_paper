# **********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : generateHeatmapCSVFgsea
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : 
#
# generateHeatmapCSVfgsea - generates an ordered heatmap based on the leadingEdge results of
# an fgsea analysis and converts them into a easy to interpret heatmap. Takes as an input an ordered csv
#from processDatafgseaCSV.R
#
# 
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      start
#
#
# **********************************************************************


generateHeatmapCSVfgsea<-function(tableGSEA,nameofplot){

library(readr)
tableGSEA<-data.matrix(tableGSEA)

library(RColorBrewer)
hmcol<-rev(brewer.pal(11,"RdBu"))

plot.name <- paste(paste(nameofplot), "plot.svg", sep="_")

plot.namePNG <- paste(paste(nameofplot), "plot.png", sep="_")



#hr <- hclust(as.dist(1-cor(t(tableGSEA), method="spearman")), method="complete")
#hc <- hclust(as.dist(1-cor(tableGSEA, method="pearson")), method="complete") 
## Tree cutting
#mycl <- cutree(hr, h=max(hr$height)/1.5); 

#mycl <- cutree(hr, k=12); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 


# mycolhc <- cm.colors(length(unique(mycl))); mycolhc <- mycolhc[as.vector(mycl)] 



#png(plot.name, height=900, width=2200)
svg(plot.name, height=8, width=25)

library(gplots)

heatmap.2(tableGSEA,
          dendrogram='none',
          col=hmcol,
          Rowv=FALSE,
          Colv=FALSE,
          trace='none',
          keysize=0.9,
          cexRow=0.8 + 1/log10(dim(tableGSEA)[1]),
          cexCol=0.3 + 1/log10(dim(tableGSEA)[1]),
          margins=c(8,40),
          srtCol=90,
          sepwidth=c(0.05,0.05)
)

dev.off()


png(plot.namePNG, height=900, width=2200)

library(gplots)

heatmap.2(tableGSEA,
          dendrogram='none',
          col=hmcol,
          Rowv=FALSE,
          Colv=FALSE,
          trace='none',
          keysize=0.9,
          cexRow=0.8 + 1/log10(dim(tableGSEA)[1]),
          cexCol=0.3 + 1/log10(dim(tableGSEA)[1]),
          margins=c(8,40),
          srtCol=90,
          sepwidth=c(0.05,0.05)
)

dev.off()
}