
#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : MainScript_3_createOrderedHeatmaps
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Creates a heatmap based upon the Normalized Enrichment Score for the Hallmark/expressionDataset
# fgsea results and Hallmark/expressionDataset_NoCutoff fgsea results. 
#
# 
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************


#imports functions
source("Functions/processDatafgseaCSV.R")
source("Functions/generateHeatmapCSVFgsea.R")

#imports packages
library(readr)



# _____________HEATMAP for  HallmarkexpressionDataset


#import fgsea table and process for heatmap 
table <-  read_csv("expressionDataset/Hallmark/HallmarkexpressionDatasettable.csv")


#process table for heatmap
tableforheatmap<-processgenesleadingedge(table)

#save table
write.table(tableforheatmap, file = "expressionDataset/Hallmark/HallmarkexpressionDatasetTableHeatmap.csv", sep = ",")


#create heatmap
output<-generateHeatmapCSVfgsea(tableforheatmap, "expressionDataset/Hallmark/HallmarkexpressionDatasetHeatmap")






# _____________HEATMAP for  HallmarkexpressionDataset_NoCutoff


#import fgsea table and process for heatmap 
table <-  read_csv("expressionDataset_NoCutoff/Hallmark/HallmarkexpressionDataset_NoCutofftable.csv")


#process table for heatmap
tableforheatmap<-processgenesleadingedge(table)

#save table
write.table(tableforheatmap, file = "expressionDataset_NoCutoff/Hallmark/HallmarkexpressionDataset_NoCutoffTableHeatmap.csv", sep = ",")


#create heatmap
output<-generateHeatmapCSVfgsea(tableforheatmap, "expressionDataset_NoCutoff/Hallmark/HallmarkexpressionDataset_NoCutoffHeatmap")



