gs_getGeneNr <- function(mod) {
  source("calc_seq_overlap.R")
  
  gcomb <- combn(mod@allGenes, m = 2)
  
  dt <- data.table(gene1 = gcomb[1,], gene2 = gcomb[2,])
  
  dt[,start1 := as.numeric(gsub(":.*$","",gsub("^.*_","", gene1)))]
  dt[,end1   := as.numeric(gsub("^.*:","",gsub("^.*_","", gene1)))]
  dt[,start2 := as.numeric(gsub(":.*$","",gsub("^.*_","", gene2)))]
  dt[,end2   := as.numeric(gsub("^.*:","",gsub("^.*_","", gene2)))]
  
  dt[,contig1 := gsub("_[0-9]*:.*$","", gene1)]
  dt[,contig2 := gsub("_[0-9]*:.*$","", gene2)]
  
  dt[,length1 := abs(start1-end1)]
  dt[,length2 := abs(start2-end2)]
  
  dt[, ol := calc_seq_overlap(start1, end1, start2, end2, contig1, contig2), by = c("gene1", "gene2")]
  
  i <- 1
  while(dt[,max(ol)] >= 0.8) {
    cat("\r",i, "gene(s) redundant by overlap.")
    dt <- dt[order(-ol)]
    
    dlg <- dt[1, gene2]
    if(dt[1,length1] < dt[1, length2])
      dlg <- dt[1, gene1]
    
    dt <- dt[gene1 != dlg & gene2 != dlg]
    i <- i + 1
  }
  cat("\n")
  
  ng <- length(unique(c(dt$gene1, dt$gene2)))
  
  return(ng)
}

#gs_getGeneNr(bsub.gs)


