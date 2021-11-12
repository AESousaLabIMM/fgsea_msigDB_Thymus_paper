
#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : MainScript_0_StartHere.R
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : The necessary libraries to run this project
#
# 
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************

#devtools
install.packages('devtools')
#tidyverse
install.packages('tidyverse')
#data.table
install.packages('data.table')
#ggplot2
install.packages('ggplot2')
#readr
install.packages('readr')

#BiocManager
install.packages('BiocManager')

#fgsea
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")


