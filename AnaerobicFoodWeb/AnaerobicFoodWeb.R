library(BacArena)
library(data.table)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
library(Biostrings)
library(stringr)
library(ggplot2)

seedEX2bigg <- function(l){
  is.ex <- str_starts(l[[1]], "EX_")
  dict <- fread("dat/SEED2VMH_translation_edited.csv", header = F)
  if( is.ex ){
    idx  <- match(gsub("_e0$","",l), gsub("\\(e\\)$","",dict$V1))  
    bigg <- gsub("_((L|D|R)\\(e\\))","__\\1", dict$V2[idx])
  }else{
    idx  <- match(l, str_extract(dict$V1, "cpd[0-9]+"))
    bigg <- str_extract(dict$V2[idx], "(?<=EX_).*?(?=\\(e\\))")
    bigg <- gsub("_((L|D|R)$)","__\\1", bigg)
  }
  return(bigg)
}
vmh2seed <- function(mod, double.underline=F){
  library(sybil)
  library(data.table)
  library(stringr)
  
  dic <- fread("dat/SEED2VMH_translation_edited.csv", col.names = c("seed","vmh"))
  ex <- findExchReact(mod)
  
  idx <- match(mod@react_id, dic$vmh)
  
  react_new <- ifelse( !is.na(idx), dic$seed[idx], react_id(mod))
  react_new <- gsub("\\(([a-z])\\)$", "_\\10", react_new)

  if( double.underline ) react_new <- gsub("_((L|D|R)\\(e\\))","__\\1", react_new)
  
  react_id(mod) <- react_new
  
  warning("Only reactions are changed!")
  return(mod)
}

models.gapseq    <- readRDS("models/anaerobic_foodweb-20201019.RDS")
models.modelseed <- readRDS("models/anaerobic_foodweb-modelseed.RDS")
models.carveme   <- readRDS("models/anaerobic_foodweb-carveme.RDS")

medium <- fread("dat/FT3.csv")
medium.carveme <- fread("dat/FT3_carveme.csv")

get_sim <- function(models, namespace.seed, tsteps=7){
  # remove some species 
  models <- models[-grep(c("GCF_000020605.1|GCF_000016525.1|GCF_000195895.1"), sapply(models, mod_desc))] # eubacterium rectale + archaeum (smithii, m.barkeri)
  models <- models[-grep(c("GCF_000024905.1|GCF_000014965.1|GCF_000013405.1"), sapply(models, mod_desc))] # syntrophic
  
  if( namespace.seed ){
    dic <- fread("~/uni/gapseq/dat/SEED2VMH_translation_edited.csv")
    mb2 <- vmh2seed(upgradeModelorg(readRDS("dat/M.barkeri_iAF692.RDS")))
    models <- append(models, mb2)
  }else{
    mb2 <- upgradeModelorg(readRDS("dat/M.barkeri_iAF692.RDS"))
    mb2.ex <- findExchReact(mb2)
    mb2@react_id[mb2.ex@react_pos] <- gsub("_((L|D|R)\\(e\\))","__\\1", mb2@react_id[mb2.ex@react_pos])
    models <- append(models, mb2)
  }
  
  arena <- Arena(n=20, m=20)
  for(i in seq_along(models)){
    mod <- models[[i]]
    if( "EX_cpd00221_e0" %in% mod@react_id ) mod <- rmReact(mod, react="EX_cpd00221_e0") # d-lactate
    bac   <- Bac(mod, setAllExInf=TRUE) 
    arena <- addOrg(arena, bac, amount = 5)
  }
  
  sub.ferm <- c("cpd00011", "cpd00029", "cpd00141", "cpd00159", "cpd00221", "cpd00211", "cpd00363", "cpd11640", "cpd01024", "cpd00036", "cpd00047")
  if( !namespace.seed) sub.ferm <- seedEX2bigg(sub.ferm)
  
  arena.medium <- paste0("EX_",medium[!compounds %in% sub.ferm, compounds],"_e0")
  if( !namespace.seed ) arena.medium <- seedEX2bigg(arena.medium)
  arena <- addSubs(arena, arena.medium, smax = medium[!compounds %in% sub.ferm, maxFlux]/100, unit = "mM")
  arena.medium.add <- c("EX_cpd00081_e0", "EX_cpd00443_e0", "EX_cpd00244_e0") # so3, 4abz, ni2 for M. barkeri
  if( !namespace.seed ) arena.medium.add <- seedEX2bigg(arena.medium.add)
  arena <- addSubs(arena, mediac = arena.medium.add, smax=1, unit="mM") 
  
  sim <- simEnv(arena, time=tsteps, sec_obj='mtf')
  
  return(sim)
}

sim.gapseq <- get_sim(models.gapseq, namespace.seed = T, tsteps=7)
sim.modelseed <- get_sim(models.modelseed, namespace.seed = T, tsteps=7)
sim.carveme <- get_sim(models.carveme, namespace.seed = F, tsteps=7)

saveRDS(sim.gapseq    , "dat/anaerobic-foodweb_bacarena2.RDS", compress = "xz")
saveRDS(sim.modelseed , "dat/anaerobic-foodweb-modelseed_bacarena2.RDS", compress = "xz")
saveRDS(sim.carveme   , "dat/anaerobic-foodweb-carveme_bacarena2.RDS", compress = "xz")
