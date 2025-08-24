# Define the function to be parallel analysis the drug.screen the score

drug_screen.score.parallel<-function(
    
    sigGenes_up,
    seurobj,
    odir,
    stim.tag,
    ctrl.tag, 
    condition,
    ncpus
    
  ){
  parallel_analysis <- function(i, sigGenes_up, seurobj, odir, stim.tag, ctrl.tag, condition) {
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(magrittr)
  
  output_dir <- file.path(odir, paste0("test"))
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if score.assignment is available
  if (!exists("score.assignment")) {
    stop("Function 'score.assignment' is not available.")
  }
  
  # The order of stim.tag should be consistent with the sigGenes
  cell.score <- score.assignment(
    seurobj = seurobj,
    odir = output_dir,
    marker.up = "yes",
    stim.tag = stim.tag,
    ctrl.tag = ctrl.tag,
    drug.name = names(sigGenes_up)[i],
    sigGenes_up = sigGenes_up[[i]],
    sigGenes_down = NULL
  )
  
  # Unify the all identification name for figure plot
  if (length(unique(cell.score$All_identity)) > 1) {
    cell.score$All_identity <- condition      # Change the original identity inconsistent into consistent value;
  }
  print(i)
  return(cell.score)
}

current_time <- Sys.time()
# Initialize snowfall and specify the number of CPU cores to use
sfInit(parallel = TRUE, cpus = ncpus)  # Change the number of cpus as needed

# Parallelize the loop
sfExport("score.assignment", "sigGenes_up", "seurobj", "odir", "stim.tag", "ctrl.tag", "condition")  #Export all variables to the workers

cell.scores <- sfSapply(1:length(sigGenes_up), parallel_analysis, simplify = FALSE,
                        sigGenes_up = sigGenes_up, 
                        seurobj = seurobj, 
                        odir = odir, 
                        stim.tag = stim.tag, 
                        ctrl.tag = ctrl.tag, 
                        condition = condition)

names(cell.scores)<- names(sigGenes_up)

# Stop snowfall cluster explicitly
  tryCatch({
    sfStop()
  }, error = function(e) {
    cat("Error occurred while stopping the cluster:", conditionMessage(e), "\n")
  })

return(cell.scores)
}



