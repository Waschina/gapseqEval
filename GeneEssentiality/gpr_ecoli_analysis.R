# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Script that compares the GPRs of                    #
# protein complexes between E. coli's curated model   #
# (iML1515) and gapseq's reconstruction.              #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

library(sybilSBML)
library(IRanges)
library(data.table)
library(stringr)

# e. coli gene loci
b_loci <- fread("essentiality.data/loci-ecol.csv")
b_loci[end < start, tmp := end]
b_loci[!is.na(tmp), end := start]
b_loci[!is.na(tmp), start := tmp]
b_loci[, tmp := NULL]
ir_B <- IRanges(start = b_loci$start, end = b_loci$end, names = b_loci$locus)

# ~ ~ ~ ~ ~ ~ ~ #
# CURATED (iML) #
# ~ ~ ~ ~ ~ ~ ~ #
iML.mod  <- readSBMLmod("models/curated/iML1515.xml")
iML.cplx <- data.table(rxn_bigg = iML.mod@react_id,
                       gpr_bigg = iML.mod@gpr)
seed_rxns <- str_match_all(iML.mod@react_attr[,1], "rxn[0-9]{5}")
names(seed_rxns) <- iML.cplx$rxn_bigg
iML.cplx <- iML.cplx[grepl("and", gpr_bigg)]
seed_rxns <- seed_rxns[iML.cplx$rxn_bigg]
iML.cplx$seed_rxn <- unlist(lapply(seed_rxns, FUN = function(x) if(length(x)==0) NA_character_ else sort(x)[1]))

# ~ ~ ~ ~ ~ ~ ~ #
# CARVEME       #
# ~ ~ ~ ~ ~ ~ ~ #
rCM.mod  <- readSBMLmod("models/carveme/ecol_protein.xml")
rCM.cplx <- data.table(rxn_cm   = rCM.mod@react_id,
                       gpr_cm   = rCM.mod@gpr)

rCM.loci <- fread("models/carveme/ecol_genes_pos.tsv")[, -2]
colnames(rCM.loci) <- c("mod.gene","start","end")
rCM.loci[end < start, tmp := end]
rCM.loci[!is.na(tmp), end := start]
rCM.loci[!is.na(tmp), start := tmp]
rCM.loci[, tmp := NULL]
ir_A <- IRanges(start = rCM.loci$start, end = rCM.loci$end, names = rCM.loci$mod.gene)

ol_est <- findOverlapPairs(ir_A, ir_B, minoverlap = 10)
mgtmp <- data.table(loci   = ol_est@second@NAMES,
                    gsterm = ol_est@first@NAMES,
                    ol     = pintersect(ol_est@first, ol_est@second)@width)
mgtmp <- mgtmp[order(gsterm, -ol)]
mgtmp <- mgtmp[!duplicated(gsterm)] # remove shorter overlaps with other genes

for(i in 1:nrow(mgtmp)) {
  qry <- mgtmp[i, gsterm]
  rpl <- mgtmp[i, loci]
  rCM.cplx[, gpr_cm := gsub(qry, rpl, gpr_cm, fixed = T)]
}
rCM.cplx[, gpr_cm := gsub("G_","", gpr_cm)]


# ~ ~ ~ ~ ~ ~ ~ #
# ModelSEED     #
# ~ ~ ~ ~ ~ ~ ~ #
rMS.mod  <- readSBMLmod("models/modelseed/ecol.sbml")
rMS.cplx <- data.table(rxn_ms   = rMS.mod@react_id,
                       gpr_ms   = rMS.mod@gpr)
rMS.cplx[, gpr_ms := gsub("(6666666\\.524383\\.peg\\.[0-9]+)([^0-9]|$)", "\\1X\\2", gpr_ms)]

rMS.loci <- prepare_MS_gene_table("models/modelseed/ecol_genes.gtbl", rMS.mod)[, -2]
colnames(rMS.loci) <- c("mod.gene","start","end")
rMS.loci[end < start, tmp := end]
rMS.loci[!is.na(tmp), end := start]
rMS.loci[!is.na(tmp), start := tmp]
rMS.loci[, tmp := NULL]
ir_A <- IRanges(start = rMS.loci$start, end = rMS.loci$end, names = rMS.loci$mod.gene)

ol_est <- findOverlapPairs(ir_A, ir_B, minoverlap = 10)
mgtmp <- data.table(loci   = ol_est@second@NAMES,
                    gsterm = ol_est@first@NAMES,
                    ol     = pintersect(ol_est@first, ol_est@second)@width)
mgtmp <- mgtmp[order(gsterm, -ol)]
mgtmp <- mgtmp[!duplicated(gsterm)] # remove shorter overlaps with other genes

for(i in 1:nrow(mgtmp)) {
  qry <- mgtmp[i, gsterm]
  qry <- paste0(qry,"X")
  rpl <- mgtmp[i, loci]
  rMS.cplx[, gpr_ms := gsub(qry, rpl, gpr_ms, fixed = T)]
}
rMS.cplx <- rMS.cplx[grepl("rxn", rxn_ms)]
rMS.cplx[, rxn_ms := gsub("_[e|p|c]0","",rxn_ms)]

# ~ ~ ~ ~ ~ ~ ~ #
# GAPSEQ        #
# ~ ~ ~ ~ ~ ~ ~ #
rGS.mod  <- readRDS("models/gapseq/ecol.RDS")
rGS.cplx <- data.table(rxn_gs   = rGS.mod@react_id,
                       gpr_gs   = rGS.mod@gpr)

rGS.loci <- data.table(modgenes = rGS.mod@allGenes)
rGS.loci[, start := as.numeric(gsub("_|\\:","",str_match(modgenes,"_[0-9]+\\:")))]
rGS.loci[, end   := as.numeric(str_match(modgenes,"[0-9]+$"))]
colnames(rGS.loci) <- c("mod.gene","start","end")
rGS.loci[end < start, tmp := end]
rGS.loci[!is.na(tmp), end := start]
rGS.loci[!is.na(tmp), start := tmp]
rGS.loci[, tmp := NULL]
ir_A <- IRanges(start = rGS.loci$start, end = rGS.loci$end, names = rGS.loci$mod.gene)

ol_est <- findOverlapPairs(ir_A, ir_B, minoverlap = 10)
mgtmp <- data.table(loci   = ol_est@second@NAMES,
                    gsterm = ol_est@first@NAMES,
                    ol     = pintersect(ol_est@first, ol_est@second)@width)
mgtmp <- mgtmp[order(gsterm, -ol)]
mgtmp <- mgtmp[!duplicated(gsterm)] # remove shorter overlaps with other genes

for(i in 1:nrow(mgtmp)) {
  qry <- mgtmp[i, gsterm]
  rpl <- mgtmp[i, loci]
  rGS.cplx[, gpr_gs := gsub(qry, rpl, gpr_gs, fixed = T)]
}
rGS.cplx <- rGS.cplx[grepl("rxn", rxn_gs)]
rGS.cplx[, rxn_gs := gsub("_[e|p|c]0","",rxn_gs)]
rGS.cplx[, gpr_gs := gsub(" | NC_000913.3_3187427:3189850","", gpr_gs, fixed = T)]

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# compare iML1515 with Recons #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
iML.cplx$test.cm <- NA
iML.cplx$test.ms <- NA
iML.cplx$test.gs <- NA

compare_boolean_expressions <- function(a,b) {
  #a <- "((b0077 | b3671 | b0507 | b0871) & (b0078))"
  #b <- "((b0077 | b0871) & (b0078))"
  a <- gsub("and", "&", a)
  a <- gsub("or", "|", a)
  b <- gsub("and", "&", b)
  b <- gsub("or", "|", b)
  
  g_words <- unlist(str_extract_all(c(a,b), "[[:alnum:]]+"))
  g_words <- unique(g_words)
  n <- length(g_words)
  
  # pffff, Friday evening coding... sorry
  for(k in 1:n) {
    #cat(" .",k,".")
    l_comb <- combn(1:n,k)
    
    for(z in 1:ncol(l_comb)) {
      b_vec <- rep(T,n)
      names(b_vec) <- g_words
      b_vec[l_comb[,z]] <- F
      
      a_eval <- with(as.list(b_vec), eval(parse(text = a)))
      b_eval <- with(as.list(b_vec), eval(parse(text = b)))

      if(a_eval != b_eval)
        return(F)
    }
  }
  return(T)
}


for(i in 1:nrow(iML.cplx)) {
  bg <- iML.cplx[i, rxn_bigg]
  sd <- iML.cplx[i, seed_rxn]
  #cat(bg,"(",sd,")")
  
  a <- iML.cplx[i, gpr_bigg]
  
  # CM
  #cat(" CM")
  b <- rCM.cplx[rxn_cm == bg, gpr_cm]
  if(length(b)==1) {
    iML.cplx[i, gpr_cm := b]
    if(b == "")
      iML.cplx[i, test.cm := F]
    if(b != "")
      iML.cplx[i, test.cm := compare_boolean_expressions(a,b)]
  } 
  
  if(!is.na(sd)) {
    # MS
    #cat(" MS")
    b <- rMS.cplx[rxn_ms == sd, gpr_ms]
    if(length(b)==1) {
      iML.cplx[i, gpr_ms := b]
      if(b == "")
        iML.cplx[i, test.ms := F]
      if(b != "")
        iML.cplx[i, test.ms := compare_boolean_expressions(a,b)]
    } 
    
    # GS
    #cat(" GS")
    b <- rGS.cplx[rxn_gs == sd, gpr_gs]
    if(length(b)==1) {
      iML.cplx[i, gpr_gs := b]
      if(b == "")
        iML.cplx[i, test.gs := F]
      if(b != "")
        iML.cplx[i, test.gs := compare_boolean_expressions(a,b)]
    } 
  }
  #cat("\n")
}

iML.cplx_small <- iML.cplx[!is.na(test.cm) & !is.na(test.ms) & !is.na(test.gs)]

prop.table(table(iML.cplx_small$test.cm));table(iML.cplx_small$test.cm)
prop.table(table(iML.cplx_small$test.ms));table(iML.cplx_small$test.ms)
prop.table(table(iML.cplx_small$test.gs));table(iML.cplx_small$test.gs)


# The correct prediction of GPR-associations for reactions catalysed by
# protein complexes remains a major challenge for the automated reconstruction
# of genome-scale metabolic models. The curated model of E. coli (iML1515, REF)
# comprises 313 reactions associated with protein complexes. Among those reactions,
# only 59 are shared with the three automated recons for E. coli using carveme,
# modelseed, and gapseq. Compared to the curated iML1515 model, 3 (6%) of the
# complex GPRs are identical in the carveme reconstruction, 6 (10%) for modelseed
# and 11 (19%) in gapseq.


# table for export supplements
out_tab <- iML.cplx[,.(biggID = rxn_bigg,
                       seedID = seed_rxn,
                       gpr_iML1515 = gpr_bigg,
                       gpr_carveme = gpr_cm,
                       test.cm,
                       gpr_modelseed = gpr_ms,
                       test.ms,
                       gpr_gapseq = gpr_gs,
                       test.gs)]
out_tab[, gpr_gapseq := gsub("\\&","and", gpr_gapseq)]
out_tab[, gpr_gapseq := gsub("\\|","or", gpr_gapseq)]
out_tab[, gpr_gapseq := gsub("\\((b[0-9]{4})\\)","\\1", gpr_gapseq)]

out_tab <- out_tab[order(is.na(gpr_carveme) + 
                           is.na(gpr_modelseed) + 
                                   is.na(gpr_gapseq))]

fwrite(out_tab, file = "Table_S_GPRs.csv")

