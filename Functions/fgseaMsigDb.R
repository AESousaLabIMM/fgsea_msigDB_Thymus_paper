

#**********************************************************************
# Project           : fgsea_msigDB_Thymus
#
# Program name      : fgseaMsigDB
#
# Author            : Susana Paço
#
# Date created      : 20211110
#
# Summary           : The function runGSEAonTest runs the fgsea algorithm on all the collections of the msigDB
# database
#
# Revision History  :
#
# Date        Author      Num    Summary
# 20211110    Susana      1      Created
# 
#
#**********************************************************************









runGSEAonTest<-function(stats, namenametransition){
  
  
  source("Functions/FunctionsForGSEA.R")
  
  library(fgsea)
  library(tidyverse)
  library(dplyr)
  
  toString(nametransition)
  
  # ------------------------Tests -----------------------------

  # -- C1 positional
  
  C1pos<- gmtPathways("MSigDb/c1.all.v7.4.symbols.gmt")
  C1postable<- createGSEA(stats,C1pos,"C1pos", nametransition)
  #namingcsv
  genetableused<-"C1pos"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  #generate csv
  #C1postable <- subset(C1postable, padj<0.05)
  tibble_with_lists_to_csv(C1postable, tablename)
  
  
  # -- Hallmark
  
  Hallmark<- gmtPathways("MSigDb/h.all.v7.4.symbols.gmt")
  Hallmarktable<- createGSEA(stats,Hallmark,"Hallmark", nametransition)
  #namingcsv
  genetableused<-"Hallmark"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  

  #generate csv
  #Hallmarktable <- subset(Hallmarktable, padj<0.05)
  tibble_with_lists_to_csv(Hallmarktable, tablename)

  # -- c2all
  
  c2all<- gmtPathways("MSigDb/c2.all.v7.4.symbols.gmt")
  c2alltable<- createGSEA(stats,c2all,"c2all", nametransition)
  #namingcsv
  genetableused<-"c2all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2alltable <- subset(c2alltable, padj<0.05)
  tibble_with_lists_to_csv(c2alltable, tablename)
  
  
  # -- c2cgp
  
  c2cgp<- gmtPathways("MSigDb/c2.cgp.v7.4.symbols.gmt")
  c2cgptable<- createGSEA(stats,c2cgp,"c2cgp", nametransition)
  #namingcsv
  genetableused<-"c2cgp"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2cgptable <- subset(c2cgptable, padj<0.05)
  tibble_with_lists_to_csv(c2cgptable, tablename)
  
  
  # -- c2can
  
  c2can<- gmtPathways("MSigDb/c2.cp.v7.4.symbols.gmt")
  c2cantable<- createGSEA(stats,c2can,"c2can", nametransition)
  #namingcsv
  genetableused<-"c2can"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2cantable <- subset(c2cantable, padj<0.05)
  tibble_with_lists_to_csv(c2cantable, tablename)
  
  # -- c2biocarta
  
  c2biocarta<- gmtPathways("MSigDb/c2.cp.biocarta.v7.4.symbols.gmt")
  c2biocartatable<- createGSEA(stats,c2biocarta,"c2biocarta", nametransition)
  #namingcsv
  genetableused<-"c2biocarta"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2biocartatable <- subset(c2biocartatable, padj<0.05)
  tibble_with_lists_to_csv(c2biocartatable, tablename)
  
  # -- c2kegg
  
  c2kegg<- gmtPathways("MSigDb/c2.cp.biocarta.v7.4.symbols.gmt")
  c2keggtable<- createGSEA(stats,c2kegg,"c2kegg", nametransition)
  #namingcsv
  genetableused<-"c2kegg"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2keggtable <- subset(c2keggtable, padj<0.05)
  tibble_with_lists_to_csv(c2keggtable, tablename)
  

  # -- c2pid
  
  c2pid<- gmtPathways("MSigDb/c2.cp.pid.v7.4.symbols.gmt")
  c2pidtable<- createGSEA(stats,c2pid,"c2pid", nametransition)
  #namingcsv
  genetableused<-"c2pid"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2pidtable <- subset(c2pidtable, padj<0.05)
  tibble_with_lists_to_csv(c2pidtable, tablename)

  # -- c2react
  
  c2react<- gmtPathways("MSigDb/c2.cp.reactome.v7.4.symbols.gmt")
  c2reacttable<- createGSEA(stats,c2react,"c2react", nametransition)
  #namingcsv
  genetableused<-"c2react"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #c2reacttable <- subset(c2reacttable, padj<0.05)
  tibble_with_lists_to_csv(c2reacttable, tablename)
  
  
  # -- C3all
  
  C3all<- gmtPathways("MSigDb/c3.all.v7.4.symbols.gmt")
  C3alltable<- createGSEA(stats,C3all,"C3all", nametransition)
  #namingcsv
  genetableused<-"C3all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3alltable <- subset(C3alltable, padj<0.05)
  tibble_with_lists_to_csv(C3alltable, tablename)
  
  
  # -- C3mir
  
  C3mir<- gmtPathways("MSigDb/c3.mir.v7.4.symbols.gmt")
  C3mirtable<- createGSEA(stats,C3mir,"C3mir", nametransition)
  #namingcsv
  genetableused<-"C3mir"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3mirtable <- subset(C3mirtable, padj<0.05)
  tibble_with_lists_to_csv(C3mirtable, tablename)
  
  
  # -- C3mirtarget
  
  C3mirtarget<- gmtPathways("MSigDb/c3.mir.mirdb.v7.4.symbols.gmt")
  C3mirtargettable<- createGSEA(stats,C3mirtarget,"C3mirtarget", nametransition)
  #namingcsv
  genetableused<-"C3mirtarget"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3mirtargettable <- subset(C3mirtargettable, padj<0.05)
  tibble_with_lists_to_csv(C3mirtargettable, tablename)
  
  
  # -- C3mirlegacy
  
  C3mirlegacy<- gmtPathways("MSigDb/c3.mir.mir_legacy.v7.4.symbols.gmt")
  C3mirlegacytable<- createGSEA(stats,C3mirlegacy,"C3mirlegacy", nametransition)
  #namingcsv
  genetableused<-"C3mirlegacy"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3mirlegacytable <- subset(C3mirlegacytable, padj<0.05)
  tibble_with_lists_to_csv(C3mirlegacytable, tablename)
  
  # -- C3tft
  
  C3tft<- gmtPathways("MSigDb/c3.tft.v7.4.symbols.gmt")
  C3tfttable<- createGSEA(stats,C3tft,"C3tft", nametransition)
  #namingcsv
  genetableused<-"C3tft"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3tfttable <- subset(C3tfttable, padj<0.05)
  tibble_with_lists_to_csv(C3tfttable, tablename)
  
  
  # -- C3tftgtrd
  
  C3tftgtrd<- gmtPathways("MSigDb/c3.tft.gtrd.v7.4.symbols.gmt")
  C3tftgtrdtable<- createGSEA(stats,C3tftgtrd,"C3tftgtrd", nametransition)
  #namingcsv
  genetableused<-"C3tftgtrd"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3tftgtrdtable <- subset(C3tftgtrdtable, padj<0.05)
  tibble_with_lists_to_csv(C3tftgtrdtable, tablename)
  
  # -- C3tftlegacy
  
  C3tftlegacy<- gmtPathways("MSigDb/c3.tft.tft_legacy.v7.4.symbols.gmt")
  C3tftlegacytable<- createGSEA(stats,C3tftlegacy,"C3tftlegacy", nametransition)
  #namingcsv
  genetableused<-"C3tftlegacy"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C3tftlegacytable <- subset(C3tftlegacytable, padj<0.05)
  tibble_with_lists_to_csv(C3tftlegacytable, tablename)
 
  # -- C4all
  
  C4all<- gmtPathways("MSigDb/c4.all.v7.4.symbols.gmt")
  C4alltable<- createGSEA(stats,C4all,"C4all", nametransition)
  #namingcsv
  genetableused<-"C4all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C4alltable <- subset(C4alltable, padj<0.05)
  tibble_with_lists_to_csv(C4alltable, tablename)
  
  # -- C4cgn
  
  C4cgn<- gmtPathways("MSigDb/c4.cgn.v7.4.symbols.gmt")
  C4cgntable<- createGSEA(stats,C4cgn,"C4cgn", nametransition)
  #namingcsv
  genetableused<-"C4cgn"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C4cgntable <- subset(C4cgntable, padj<0.05)
  tibble_with_lists_to_csv(C4cgntable, tablename)
  
  # -- C4cm
  
  C4cm<- gmtPathways("MSigDb/c4.cm.v7.4.symbols.gmt")
  C4cmtable<- createGSEA(stats,C4cm,"C4cm", nametransition)
  #namingcsv
  genetableused<-"C4cm"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C4cmtable <- subset(C4cmtable, padj<0.05)
  tibble_with_lists_to_csv(C4cmtable, tablename)
  
  # -- C5all
  
  C5all<- gmtPathways("MSigDb/c5.all.v7.4.symbols.gmt")
  C5alltable<- createGSEA(stats,C5all,"C5all", nametransition)
  #namingcsv
  genetableused<-"C5all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C5alltable <- subset(C5alltable, padj<0.05)
  tibble_with_lists_to_csv(C5alltable, tablename)
  
  # -- C5bp
  
  C5bp<- gmtPathways("MSigDb/c5.go.bp.v7.4.symbols.gmt")
  C5bptable<- createGSEA(stats,C5bp,"C5bp", nametransition)
  #namingcsv
  genetableused<-"C5bp"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C5bptable <- subset(C5bptable, padj<0.05)
  tibble_with_lists_to_csv(C5bptable, tablename)
  
  # -- C5cc
  
  C5cc<- gmtPathways("MSigDb/c5.go.cc.v7.4.symbols.gmt")
  C5cctable<- createGSEA(stats,C5cc,"C5cc", nametransition)
  #namingcsv
  genetableused<-"C5cc"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C5cctable <- subset(C5cctable, padj<0.05)
  tibble_with_lists_to_csv(C5cctable, tablename)
  
  # -- C5mf
  
  C5mf<- gmtPathways("MSigDb/c5.go.mf.v7.4.symbols.gmt")
  C5mftable<- createGSEA(stats,C5mf,"C5mf", nametransition)
  #namingcsv
  genetableused<-"C5mf"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C5mftable <- subset(C5mftable, padj<0.05)
  tibble_with_lists_to_csv(C5mftable, tablename)
  
  # -- C6all
  
  C6all<- gmtPathways("MSigDb/c6.all.v7.4.symbols.gmt")
  C6alltable<- createGSEA(stats,C6all,"C6all", nametransition)
  #namingcsv
  genetableused<-"C6all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C6alltable <- subset(C6alltable, padj<0.05)
  tibble_with_lists_to_csv(C6alltable, tablename)
  
  
  # -- C7all
  
  C7all<- gmtPathways("MSigDb/c7.all.v7.4.symbols.gmt")
  C7alltable<- createGSEA(stats,C7all,"C7all", nametransition)
  #namingcsv
  genetableused<-"C7all"
  #name csv
  tablename<-paste(nametransition,"/",genetableused,"/",genetableused,nametransition,"table.csv", sep="")
  
  #generate csv
  #C7alltable <- subset(C7alltable, padj<0.05)
  tibble_with_lists_to_csv(C7alltable, tablename)
  
  
}