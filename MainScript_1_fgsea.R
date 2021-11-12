


#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : 1_MainScript_fgsea
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Script that runs the FGSEA algorithm on  expressionDataset.txt and 
# expressionDataset_NoCutoff.txt and outputs a csv with results for each run, as well as a sticks plot for the top
# pathways and an enrichment plot for the significant pathways
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
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
source("Functions/fgseaMsigDb.R")
source("Functions/GenerateFolderTree.R")




# _________________________________ RUN FGSEA on  expressionDataset_NoCutoff

#substitute the inout fot the appropriate table. It's important that the IDs are HGNC


# --------- generate folder structure
basedir<-getwd()

#give this the name of the input file
nametransition<-"expressionDataset_NoCutoff"

#generates the folder strcture for the output
generatefolderstructure(nametransition)

setwd(basedir)


# -----------get and prepare test data

expressionDataset_NoCutoff <- read_table("expressionDataset_NoCutoff.txt")
expressionDataset_NoCutoff<-data.frame(expressionDataset_NoCutoff$ensembl_gene_id, expressionDataset_NoCutoff$hgnc_symbol, expressionDataset_NoCutoff$logFC)
names(expressionDataset_NoCutoff)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task - genes and logFC
expressionDataset_NoCutoff<-data.frame(expressionDataset_NoCutoff$hgnc_symbol,expressionDataset_NoCutoff$logFC)
names(expressionDataset_NoCutoff)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
expressionDataset_NoCutoff<-arrange(expressionDataset_NoCutoff, desc(expressionDataset_NoCutoff$logFC))
expressionDataset_NoCutoff <- deframe(expressionDataset_NoCutoff)


runGSEAonTest(expressionDataset_NoCutoff, nametransition)



# _________________________________ RUN FGSEA on  expressionDataset

#substitute the inout fot the appropriate table. It's important that the IDs are HGNC


# --------- generate folder structure
basedir<-getwd()

#give this the name of the input file
nametransition<-"expressionDataset"

#generates the folder strcture for the output
generatefolderstructure(nametransition)

setwd(basedir)


# -----------get and prepare test data

expressionDataset<- read_table2("Data/expressionDataset.txt")
expressionDataset<-data.frame(expressionDataset$ensembl_gene_id, expressionDataset$hgnc_symbol, expressionDataset$logFC)
names(expressionDataset)<-c("ensembl_gene_id", "hgnc_symbol", "logFC")

#get just the data for the task - genes and logFC
expressionDataset<-data.frame(expressionDataset$hgnc_symbol,expressionDataset$logFC)
names(expressionDataset)<-c("hgnc_symbol", "logFC")

#arrange data for fgsea We need a named list with the fold change values arrange from highest
#to lowest
expressionDataset<-arrange(expressionDataset, desc(expressionDataset$logFC))
expressionDataset <- deframe(expressionDataset)


runGSEAonTest(expressionDataset, nametransition)
