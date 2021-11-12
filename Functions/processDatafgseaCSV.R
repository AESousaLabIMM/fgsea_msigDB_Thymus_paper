# **********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : processDatafgseaCSV
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Process csv from fgsea into usable data for
# the heatmap
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      start
#
#
# **********************************************************************





processgenesleadingedge<-function(tableGSEA){

library(tidyverse)  
#remove unnecessary tables
#tableGSEA$X1<-NULL
  

tableGSEA$pval<-NULL
tableGSEA$ES<-NULL
tableGSEA$nMoreExtreme<-NULL
tableGSEA$size<-NULL


tableGSEA$pathway<-tableGSEA$pathway%>% str_replace("HALLMARK_", "")

#filter for padj<0.05
#tableGSEA<-tableGSEA %>% filter(padj < 0.05)



tableGSEA$padj<-NULL

tableGSEA <- tableGSEA[order(-tableGSEA$NES),]

#convert leadingEdge into dummy variables

tableGSEA$leadingEdge <- as.character(tableGSEA$leadingEdge)
resp.split <- strsplit(tableGSEA$leadingEdge, ",")

lev <- unique(unlist(resp.split))

resp.dummy <- lapply(resp.split, function(x) table(factor(x, levels=lev)))

tableGSEA2 <- with(tableGSEA, data.frame(pathway, NES, do.call(rbind, resp.dummy)))
tableGSEA2
tableGSEA<-tableGSEA2


# set pathway as index
rownames(tableGSEA) <- tableGSEA$pathway
tableGSEA$pathway<- NULL


#extract NES
NESvalues<-tableGSEA$NES
tableGSEA$NES<-NULL
#multiply NES for the dummy variables 
tableGSEA<-tableGSEA*NESvalues



return(tableGSEA)
}

