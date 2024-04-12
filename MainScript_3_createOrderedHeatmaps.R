
#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : MainScript_3_createOrderedHeatmaps
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : Creates a heatmap based upon the Normalized Enrichment Score for the Hallmark/TregvsTconv_Thy_DEG
# fgsea results and Hallmark/TregvsTconv_Thy_DEGnoco fgsea results. 
#
# 
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 20240412    Susana      1      Updated
# 
#
#**********************************************************************


#imports functions
source("Functions/processDatafgseaCSV.R")
source("Functions/generateHeatmapCSVFgsea.R")

#imports packages
library(readr)



# _____________HEATMAP for  HallmarkTregvsTconv_Thy_DEG


#import fgsea table and process for heatmap 
table <-  read_csv("TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEGtable.csv")


#process table for heatmap
tableforheatmap<-processgenesleadingedge(table)

#save table
write.table(tableforheatmap, file = "TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEGTableHeatmap.csv", sep = ",")


#create heatmap
output<-generateHeatmapCSVfgsea(tableforheatmap, "TregvsTconv_Thy_DEG/Hallmark/HallmarkTregvsTconv_Thy_DEGHeatmap")






# _____________HEATMAP for  HallmarkTregvsTconv_Thy_DEGnoco


#import fgsea table and process for heatmap 
table <-  read_csv("TregvsTconv_Thy_DEGnoco/Hallmark/HallmarkTregvsTconv_Thy_DEGnocotable.csv")


#process table for heatmap
tableforheatmap<-processgenesleadingedge(table)

#save table
write.table(tableforheatmap, file = "TregvsTconv_Thy_DEGnoco/Hallmark/HallmarkTregvsTconv_Thy_DEGnocoTableHeatmap.csv", sep = ",")


#create heatmap
output<-generateHeatmapCSVfgsea(tableforheatmap, "TregvsTconv_Thy_DEGnoco/Hallmark/HallmarkTregvsTconv_Thy_DEGnocoHeatmap")



