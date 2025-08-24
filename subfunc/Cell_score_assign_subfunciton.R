score.assignment<- function(
    seurobj,
    odir,
    ctrl.tag,
    stim.tag,
    marker.up,
    drug.name,
    sigGenes_up,
    sigGenes_down=NULL
){
  
  if(marker.up=="yes"){
    
    sigGenes_filtered<-sigGenes_up[sigGenes_up %in% row.names(seurobj)]
    gene_expression <- GetAssayData(object = seurobj@assays$RNA, layer = "counts")[sigGenes_filtered, ]
    
    # Count the number of up-marker genes(UMI>0)
    #num_genes_umi_gt_0 <- colSums(gene_expression > 0)
    
    ## count the sum UMIs of up-marker genes
    sum_umi_genes <-  colSums(gene_expression)
    seurobj[["sigGene_sum_umi"]] <- sum_umi_genes
    #seurobj[["sigGene_num"]]<-num_genes_umi_gt_0
    SigGene_UMI<-seurobj[[]]
    
    SigGene_UMI<-SigGene_UMI %>% rename("All_identity"="orig.ident","All_Count_RNA"="nCount_RNA","All_Feature_RNA"="nFeature_RNA")
    
    SigGene_UMI[,"UMI_ratio"]<-SigGene_UMI[,"sigGene_sum_umi"]/SigGene_UMI[,"All_Count_RNA"]
    #SigGene_UMI[,"NUM_ratio"]<-SigGene_UMI[,"sigGene_num"]/length(sigGenes_up)
    #SigGene_UMI[,"cell_score"]<-SigGene_UMI[,"UMI_ratio"]*SigGene_UMI[,"NUM_ratio"]*1000
    
    SigGene_UMI[,"Barcode"]<-row.names(SigGene_UMI)
    SigGene_UMI<-SigGene_UMI[,c("Barcode","All_identity", "Best_Tag", "UMI_ratio")]
    
    #head(SigGene_UMI);dim(SigGene_UMI)
    
  }else{
    
  #  seurobj<-subset(seurobj, nCount_RNA>=min.num.umi & nFeature_RNA>= min.num.gene & nFeature_RNA<=max.num.gene & Ratio_Best_Tag_In_Drop>=BT.ratio & percent.mt<=pct.mt)
    
    sigGenes_filtered_up<-sigGenes_up[sigGenes_up %in% row.names(seurobj)]
    sigGenes_filtered_down<-sigGenes_down[sigGenes_down %in% row.names(seurobj)]
    gene_expression_up <- GetAssayData(object = seurobj@assays$RNA, layer = "counts")[sigGenes_filtered_up, ]
    gene_expression_down <- GetAssayData(object = seurobj@assays$RNA, layer = "counts")[sigGenes_filtered_down, ]
    sum_umi_genes_up <-  colSums(gene_expression_up)
    sum_umi_genes_down <-  colSums(gene_expression_down)
    seurobj[["sigGene_sum_umi_up"]] <- sum_umi_genes_up
    seurobj[["sigGene_sum_umi_down"]] <- sum_umi_genes_down
    seurobj[["sigGene_sum_umi"]] <- seurobj[["sigGene_sum_umi_up"]] - seurobj[["sigGene_sum_umi_down"]] 
    
    SigGene_UMI<-seurobj[[]]
    
    SigGene_UMI<-SigGene_UMI %>% rename("All_identity"="orig.ident","All_Count_RNA"="nCount_RNA","All_Feature_RNA"="nFeature_RNA")
    
    SigGene_UMI[,"UMI_ratio"]<-SigGene_UMI[,"sigGene_sum_umi"]/SigGene_UMI[,"All_Count_RNA"]
    SigGene_UMI[,"Barcode"]<-row.names(SigGene_UMI)
    
    #head(SigGene_UMI);dim(SigGene_UMI)
  }
  
  # violin plot all cells treatment vs control 
  ggplot(SigGene_UMI, aes(x=All_identity, y=UMI_ratio)) +
    theme_classic() + 
    geom_violin(data = SigGene_UMI[SigGene_UMI$Best_Tag==ctrl.tag,], trim=TRUE, fill=NA, aes(color = "grey"))+
    geom_violin(data = SigGene_UMI[SigGene_UMI$Best_Tag==stim.tag,], trim=TRUE, fill=NA, aes(color = "red"))+
    # scale_x_discrete(limits=c("Antibody"))+
    labs(x="Drug", y="CellSign", title="") + 
    scale_colour_manual(name="Treatment", 
                        labels=c(ctrl.tag, stim.tag), 
                        values=c("grey"="grey","red"="red"))
  ggsave(filename = file.path(odir, paste0(drug.name,"_violin_all_cell_CellSign.png")) ,bg='transparent', width = 8, height = 6)
  
  # AUC curve
  
  SigGene_UMI[,"Label"] <- 1; SigGene_UMI[SigGene_UMI$Best_Tag==ctrl.tag,"Label"] <- 0
  
  par(pty="s")
  pdf(paste0(odir,"/",drug.name,"_ROC_AUC.pdf",sep = ""), height = 5, width = 6)
  
  myroc<-roc(SigGene_UMI$Label, SigGene_UMI$UMI_ratio, direction="<", plot=TRUE, percent = TRUE, legacy.axes=TRUE, 
             xlab="False Positive Precentage", ylab = "True Positive Percentage",
             col= "red", lwd =2,
             #   print.auc = TRUE, print.auc.x=95, print.auc.y=100
  )
  
  legend("bottomright", legend=c(paste0("AUC: ",round(myroc$auc[[1]],2), "%")),col= c("red"), lwd = 2)
  dev.off()
  
  SigGene_UMI[,"AUC"]<-round(myroc$auc,2)
  
  return(SigGene_UMI) 
}
