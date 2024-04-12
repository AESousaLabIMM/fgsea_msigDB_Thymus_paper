

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
# in the fgsea algorithm run in TregvsTconv_Thy_DEG.txt with the Hallmark collection of MsigDB
# 
# The results get stored in the same folder as the other Hallmark/TregvsTconv_Thy_DEG.txt fgsea run
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 20240412    Susana      1      Updated
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

TregvsTconv_Thy_DEG<- read_table2("Data/TregvsTconv_Thy_DEG.txt")
TregvsTconv_Thy_DEG<-data.frame(TregvsTconv_Thy_DEG$ensembl_gene_id, TregvsTconv_Thy_DEG$hgnc_symbol, TregvsTconv_Thy_DEG$logFC)
names(TregvsTconv_Thy_DEG)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task
TregvsTconv_Thy_DEG<-data.frame(TregvsTconv_Thy_DEG$hgnc_symbol,TregvsTconv_Thy_DEG$logFC)
names(TregvsTconv_Thy_DEG)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
TregvsTconv_Thy_DEG<-arrange(TregvsTconv_Thy_DEG, desc(TregvsTconv_Thy_DEG$logFC))
TregvsTconv_Thy_DEG <- deframe(TregvsTconv_Thy_DEG)




#plots

HALLMARK_IL2_STAT5_SIGNALINGname<-"TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEG_HALLMARK_IL2_STAT5_SIGNALING.pdf"

pdf(file = HALLMARK_IL2_STAT5_SIGNALINGname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
               TregvsTconv_Thy_DEG, ticksSize = 0.8) + labs(title="IL2 STAT5 SIGNALING")
dev.off()


HALLMARK_TNFA_SIGNALING_VIA_NFKBname<-"TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEG_HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf"

pdf(file = HALLMARK_TNFA_SIGNALING_VIA_NFKBname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               TregvsTconv_Thy_DEG, ticksSize = 0.8) + labs(title="TNFA SIGNALING VIA NFKB")
dev.off()

HALLMARK_INFLAMMATORY_RESPONSEname<-"TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEG_HALLMARK_INFLAMMATORY_RESPONSE.pdf"

pdf(file = HALLMARK_INFLAMMATORY_RESPONSEname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
               TregvsTconv_Thy_DEG, ticksSize = 0.8) + labs(title="INFLAMMATORY RESPONSE")
dev.off()


HALLMARK_GLYCOLYSISname<-"TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEG_HALLMARK_GLYCOLYSIS.pdf"

pdf(file = HALLMARK_GLYCOLYSISname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_GLYCOLYSIS"]],
               TregvsTconv_Thy_DEG, ticksSize = 0.8) + labs(title="GLYCOLYSIS")
dev.off()


HALLMARK_ESTROGEN_RESPONSE_LATEname<-"TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEG_HALLMARK_ESTROGEN_RESPONSE_LATE.pdf"

pdf(file = HALLMARK_ESTROGEN_RESPONSE_LATEname,h=10,w=12)
plotEnrichmentfgsea(Hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],
               TregvsTconv_Thy_DEG, ticksSize = 0.8) + labs(title="ESTROGEN RESPONSE LATE")
dev.off()
