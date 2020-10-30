prediction.stats <- function(conf_mat) {
  out <- c()
  
  if(ncol(conf_mat)==1)
    conf_mat <- cbind(matrix(c(0,0),ncol = 1),conf_mat)
  
  out['F1.score']    <- 2*conf_mat[2,2] / (2*conf_mat[2,2] + conf_mat[2,1] + conf_mat[1,2])
  out['Precision']   <- conf_mat[2,2] / (conf_mat[2,2] + conf_mat[1,2])
  out['Sensitivity'] <- conf_mat[2,2] / (conf_mat[2,2] + conf_mat[2,1])
  out['Specificity'] <- conf_mat[1,1] / (conf_mat[1,1] + conf_mat[1,2])
  out['Accuracy']    <- (conf_mat[1,1] +  conf_mat[2,2]) / (conf_mat[1,1] + conf_mat[1,2] + conf_mat[2,1] + conf_mat[2,2])
    
  return(out)
}
