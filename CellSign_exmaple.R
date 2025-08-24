library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(caret)
library(pROC)
library(snowfall)
library(magrittr)

project_dir<- getwd()
script_dir<-file.path(project_dir, "subfunc/")
odir<-file.path(project_dir, "Figs/")
if(!dir.exists(odir)) { dir.create(odir, recursive = T) }
Data_dir<-file.path(project_dir, "Dataset/")

source(paste0(script_dir, "Cell_score_assign_subfunciton.R"))
source(paste0(script_dir, "parallel_drug_screen_score.R"))
source(paste0(script_dir, "drug_screen_plot.R"))

# read the single cell expression matrix: 10X data
MCF_B1<- readRDS(paste0(Data_dir, "MCF7_HF2FWBGXC.rds"))

seurobj.sub<-subset(MCF_B1, nCount_RNA>=1000 & nFeature_RNA >= 1000 & nFeature_RNA<=10000 & Ratio_Best_Tag_In_Drop >=0.9 & percent.mt <=15)

linc.marker.list<-readRDS(paste0(Data_dir, "linc.MCF7.niclosamide.individual.RDS"))

## 1) calculate MCF7 B1 niclosamide treatment
stim.tag<-c("MCF7-Niclosamide"); ctrl.tag<-"MCF7-DMF"; condition="MCF7_B1"

drug_screen.screen <- drug_screen.plot(
  sigGenes_up=linc.marker.list,
  seurobj=seurobj.sub,
  odir=file.path(odir, condition, "All_drugs"),
  stim.tag=stim.tag,
  ctrl.tag=ctrl.tag, 
  condition=condition,
  ncpus=7
)

saveRDS(drug_screen.screen, paste0(odir, "all_drug_cell_scores_MCF7_B1.rds"))


