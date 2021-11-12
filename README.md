# fgsea_msigDB_Thymus_paper



The script should be run as follows:

1 - Run MainScript_0_StartHere.R to install necessary packages
2 - Run MainScript_1_fgsea.R to run the fgsea. 
  it runs in the expression dataset here named "expressionDataset" (replace the expression by the name of the file) and the expression dataset without previous cutoff here named "expressionDataset_NoCutoff" (replace the expression by the name of the file). The results are stored in a folder tree generated automatically we the name of the file you give it as root.
3 - Run   MainScript_3_createOrderedHeatmaps.R to convert the csv outputs of the Hallmark/expression dataset fgsea into a heatmap. The results will be stored in the hallmark folder in its respective results folder. 