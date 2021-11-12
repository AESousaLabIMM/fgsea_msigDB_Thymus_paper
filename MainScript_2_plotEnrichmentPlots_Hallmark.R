

#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : MainScipt_2_plotEnrichmentPlots_Hallmark.R

#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : scripts that create the enrichment plots for the significant pathways
# in the fgsea algorithm run in expressionDataset.txt with the Hallmark collection of MsigDB
# 
# The results get stored in the same folder as the other Hallmark/expressionDataset.txt fgsea run
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************








# retrieve hallmakr enrichments



# imports

library(devtools)
library(tidyverse)
library(data.table)
library(fgsea)
library(ggplot2)
library(readr)

#get functions to perform GSEA
source("Functions/FunctionsForGSEA.R")
source("Functions/plotEnrichmentfgsea.R")


#pathways
Hallmark<- gmtPathways("MsigDB/h.all.v7.4.symbols.gmt")


#stats

# -----------get and prepare test data

expressionDataset<- read_table2("Data/expressionDataset.txt")
expressionDataset<-data.frame(expressionDataset$ensembl_gene_id, expressionDataset$hgnc_symbol, expressionDataset$logFC)
names(expressionDataset)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task
expressionDataset<-data.frame(expressionDataset$hgnc_symbol,expressionDataset$logFC)
names(expressionDataset)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
expressionDataset<-arrange(expressionDataset, desc(expressionDataset$logFC))
expressionDataset <- deframe(expressionDataset)




#plots

HALLMARK_IL2_STAT5_SIGNALINGname<-"expressionDataset/Hallmark/HallmarkexpressionDataset_HALLMARK_IL2_STAT5_SIGNALING.pdf"

pdf(file = HALLMARK_IL2_STAT5_SIGNALINGname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
               expressionDataset, ticksSize = 0.8) + labs(title="IL2 STAT5 SIGNALING")
dev.off()


HALLMARK_TNFA_SIGNALING_VIA_NFKBname<-"expressionDataset/Hallmark/HallmarkexpressionDataset_HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf"

pdf(file = HALLMARK_TNFA_SIGNALING_VIA_NFKBname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               expressionDataset, ticksSize = 0.8) + labs(title="TNFA SIGNALING VIA NFKB")
dev.off()

HALLMARK_INFLAMMATORY_RESPONSEname<-"expressionDataset/Hallmark/HallmarkexpressionDataset_HALLMARK_INFLAMMATORY_RESPONSE.pdf"

pdf(file = HALLMARK_INFLAMMATORY_RESPONSEname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
               expressionDataset, ticksSize = 0.8) + labs(title="INFLAMMATORY RESPONSE")
dev.off()


HALLMARK_GLYCOLYSISname<-"expressionDataset/Hallmark/HallmarkexpressionDataset_HALLMARK_GLYCOLYSIS.pdf"

pdf(file = HALLMARK_GLYCOLYSISname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_GLYCOLYSIS"]],
               expressionDataset, ticksSize = 0.8) + labs(title="GLYCOLYSIS")
dev.off()


HALLMARK_ESTROGEN_RESPONSE_LATEname<-"expressionDataset/Hallmark/HallmarkexpressionDataset_HALLMARK_ESTROGEN_RESPONSE_LATE.pdf"

pdf(file = HALLMARK_ESTROGEN_RESPONSE_LATEname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],
               expressionDataset, ticksSize = 0.8) + labs(title="ESTROGEN RESPONSE LATE")
dev.off()
