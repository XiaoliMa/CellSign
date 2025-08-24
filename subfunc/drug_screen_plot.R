drug_screen.plot <- function(
  
  sigGenes_up,
  seurobj,
  odir,
  stim.tag,
  ctrl.tag, 
  condition,
  ncpus=ncpus

  ){

    drug.screen.scores <- drug_screen.score.parallel(
        sigGenes_up=sigGenes_up,
        seurobj=seurobj,
        odir=odir,
        stim.tag=stim.tag,
        ctrl.tag=ctrl.tag, 
        condition=condition,
        ncpus=ncpus
  )
    
    # plot one drug treatment
    if(length(unique(drug.screen.scores[[1]]$Best_Tag))<3){

    name.auc<-data.frame()
    for(i in c(1:length(drug.screen.scores))){
      name.auc.tmp <- data.frame(
        drug.name = names(drug.screen.scores)[i],
        auc.value = drug.screen.scores[[i]]$AUC[1],
        drug=stim.tag
      )
      name.auc<- rbind(name.auc,name.auc.tmp)
    }
    
  name.auc<-name.auc[order(name.auc$auc.value, decreasing = T),]
  top5 <- name.auc %>% group_by(drug) %>% top_n(5, auc.value)
  non.top5.auc <- anti_join(name.auc, top5, by = c("drug", "auc.value"))
  data_clean <- na.omit(name.auc) 
  
  p<-ggplot(name.auc, aes(x = drug, y = auc.value)) + theme_minimal(base_size = 12) + theme_classic() +
    stat_boxplot(geom = "errorbar",width = 0.15) + 
    geom_boxplot(outliers = FALSE) +
    geom_point(position = position_jitter(width = 0.25, height = 0.01), data = non.top5.auc, color = "grey") +
    geom_point(data = top5, aes(x = drug, y = auc.value), color = c("red"), size =2) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # Center and bold the title
    ggtitle("Drug Screening") + ylim(c(min(name.auc$auc.value)-0.5,100))+
    ylab("S(cs, up)") + xlab("Drug") 
  
  for(i in seq_along(top5$drug.name)){
    
    y.position <- top5$auc.value[i] - 3 * i
    x.position <- c(1, 0.8, 1.2, 1, 0.8)
    p <- p+geom_text(label = top5$drug.name[i], x = x.position[i], y= y.position , vjust = -1.5, color = "red") # Add text label for median
  }
  # Add all text layers to the plot
  p 
  ggsave(filename = file.path(odir, paste0(condition,"/all_drug_screen_boxplot_AUC_CellSign.png")) ,bg='transparent', width = 8, height = 6, dpi = 300)
  return(list("drug.screen.scores"=drug.screen.scores,"name.auc"=name.auc))
  
    }
  
}
