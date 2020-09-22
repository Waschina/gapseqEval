predGeneEssentiality.GS <- function(mod, loci.file, gess.file) {
  multi.contig.mode <- F
  loci <- fread(loci.file)
  loci[, gs.kogr := NA_real_]
  
  gess <- fread(gess.file)
  loci <- loci[locus %in% gess$gene]
  if("genomic.element" %in% colnames(loci)) {
    multi.contig.mode <- T
    print("multi contig mode")
  }
  
  model.genes <- data.table(mod.gene = mod@allGenes)
  
  # parsing gapseq-style "gene names"
  model.genes[,contig := str_replace_all(mod.gene, "_[0-9]*:[0-9]*$","")]
  model.genes[,range  := str_replace_all(mod.gene, "^.*_","")]
  model.genes[,start  := str_replace_all(range, ":[0-9]*$","")]
  model.genes[,end    := str_replace_all(range, "^[0-9]*:","")]
  
  for(i in 1:nrow(loci)) {
    cat("\r",i,"/",nrow(loci))
    if(multi.contig.mode)
      model.genes$s.ol <- calc_seq_overlap(loci[i,start],loci[i,end], model.genes$start, model.genes$end, loci[i, genomic.element], model.genes$contig)
    else
      model.genes$s.ol <- calc_seq_overlap(loci[i,start],loci[i,end], model.genes$start, model.genes$end)
    #print(model.genes[s.ol > 0.5])
    loci[i, mod.genes := paste(model.genes[s.ol > 0.25, mod.gene], collapse = ";")]
    
    n.ko <- nrow(model.genes[s.ol > 0.25])
    
    if(n.ko > 0) {
      ko.sol <- optimizeProb(mod, gene = model.genes[s.ol > 0.25, mod.gene],
                             lb = rep(0, n.ko), ub = rep(0, n.ko))
      
      loci[i, gs.kogr := ko.sol@lp_obj]
    }
    
  }
  cat("\n")
  loci[, ess.gapseq := ifelse(gs.kogr < 0.01,"yes","no")]
  
  
  # merge with currected model prediction and experimentally-determined essentiality
  loci <- merge(loci, gess, by.x = "locus", by.y = "gene")
  
  
  return(list(loci = loci,
              model.genes = model.genes))
}
