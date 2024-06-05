


#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : 1_MainScript_fgsea
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Script that runs the FGSEA algorithm on  Data/TregvsTconv_Thy_DEG.txt and 
# Data/TregvsTconv_Thy_DEGnoco.txt and outputs a csv with results for each run, as well as a sticks plot for the top
# pathways and an enrichment plot for the significant pathways
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 20240412    Susana      1      Updated
# 
#
#**********************************************************************





#necessary imports
library(devtools)
library(tidyverse)
library(data.table)
library(fgsea)
library(ggplot2)
library(readr)

#get functions to perform GSEA
source("Functions/FunctionsForGSEA.R")

#v7.4
source("Functions/fgseaMsigDb.R")



source("Functions/GenerateFolderTree.R")




# _________________________________ RUN FGSEA on  TregvsTconv_Thy_DEGnoco


# --------- generate folder structure
basedir<-getwd()

#give this the name of the input file
nametransition<-"TregvsTconv_Thy_DEGnoco"

#generates the folder strcture for the output
generatefolderstructure(nametransition)

setwd(basedir)


# -----------get and prepare test data

TregvsTconv_Thy_DEGnoco <- read_table("Data/TregvsTconv_Thy_DEGnoco.txt")
TregvsTconv_Thy_DEGnoco<-data.frame(TregvsTconv_Thy_DEGnoco$ensembl_gene_id, TregvsTconv_Thy_DEGnoco$hgnc_symbol, TregvsTconv_Thy_DEGnoco$logFC)
names(TregvsTconv_Thy_DEGnoco)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task - genes and logFC
TregvsTconv_Thy_DEGnoco<-data.frame(TregvsTconv_Thy_DEGnoco$hgnc_symbol,TregvsTconv_Thy_DEGnoco$logFC)
names(TregvsTconv_Thy_DEGnoco)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
TregvsTconv_Thy_DEGnoco<-arrange(TregvsTconv_Thy_DEGnoco, desc(TregvsTconv_Thy_DEGnoco$logFC))
TregvsTconv_Thy_DEGnoco <- deframe(TregvsTconv_Thy_DEGnoco)


runGSEAonTest(TregvsTconv_Thy_DEGnoco, nametransition)



# _________________________________ RUN FGSEA on  TregvsTconv_Thy_DEG


# --------- generate folder structure
basedir<-getwd()

#give this the name of the input file
nametransition<-"TregvsTconv_Thy_DEG"

#generates the folder strcture for the output
generatefolderstructure(nametransition)

setwd(basedir)


# -----------get and prepare test data

TregvsTconv_Thy_DEG<- read_table2("Data/TregvsTconv_Thy_DEG.txt")
TregvsTconv_Thy_DEG<-data.frame(TregvsTconv_Thy_DEG$ensembl_gene_id, TregvsTconv_Thy_DEG$hgnc_symbol, TregvsTconv_Thy_DEG$logFC)
names(TregvsTconv_Thy_DEG)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task - genes and logFC
TregvsTconv_Thy_DEG<-data.frame(TregvsTconv_Thy_DEG$hgnc_symbol,TregvsTconv_Thy_DEG$logFC)
names(TregvsTconv_Thy_DEG)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
TregvsTconv_Thy_DEG<-arrange(TregvsTconv_Thy_DEG, desc(TregvsTconv_Thy_DEG$logFC))
TregvsTconv_Thy_DEG <- deframe(TregvsTconv_Thy_DEG)


runGSEAonTest(TregvsTconv_Thy_DEG, nametransition)
