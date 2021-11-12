

#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : GenerateolderTree
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : 
#
# generatefolderstructure creates a folder tree to store the results of 
# the fgsea analysis with the msigDB database. The root folder adopts the name
# of the original data input
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************



generatefolderstructure<-function(testname){
  
  basepath<-paste(testname) 
  toString(basepath)
  dir.create(basepath)
  
  #list folders
  msigtables<- list(      "C1pos",
                          "Hallmark",
                          "c2all",
                          "c2cgp",
                          "c2can",
                          "c2biocarta",
                          "c2kegg",
                          "c2pid",
                          "c2react",
                          "C3all",
                          "C3mir",
                          "C3mirtarget",
                          "C3mirlegacy",
                          "C3tft",
                          "C3tftgtrd",
                          "C3tftlegacy",
                          "C4all",
                          "C4cgn",
                          "C4cm",
                          "C5all",
                          "C5bp",
                          "C5cc",
                          "C5mf",
                          "C6all",
                          "C7all")
  setwd(basepath)  
  
  lapply(msigtables,dir.create,recursive = TRUE)
  
}