library(data.table)
library(Biostrings)
library(stringr)
library(ggplot2)

seed <- fread("dat/seed_metabolites_edited.tsv")
genome.desc <- fread("dat/genome_desc.csv")

sim.gapseq    <- readRDS("dat/anaerobic-foodweb_bacarena2.RDS")
sim.modelseed <- readRDS("dat/anaerobic-foodweb-modelseed_bacarena.RDS")
sim.carveme   <- readRDS("dat/anaerobic-foodweb-carveme_bacarena.RDS")

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

sub.ferm.seed <- c("cpd00011", "cpd00029", "cpd00141", "cpd00159", "cpd00221", "cpd00211", "cpd00363", "cpd11640", "cpd01024", "cpd00036", "cpd00047")
sub.ferm.bigg <- seedEX2bigg(sub.ferm.seed)

sim.activity <- function(sim, namespace.seed, sub.ferm, agora=F){
  if( namespace.seed ){
    dat.sub <- unique(c(setdiff(sub.ferm, "cpd00011"), "cpd00239"))
    dat <- data.table(plotSpecActivity(sim, subs=paste0("EX_",dat.sub,"_e0"),useNames = T, rm_unused = T, ret_data = T))
    dat[,id:=seed$id[match(tolower(gsub("-e0$","",sub)), tolower(seed$name))]]
    dat[sub=="raffinose", id:="cpd00382"]
    dat[,spec.name:=trimws(genome.desc$name[match(gsub(".fna.sbml","", spec), gsub(".fna.gz","",genome.desc$file))])]
  }else{
    dat.sub <- unique(c(setdiff(sub.ferm, seedEX2bigg("cpd00011")), seedEX2bigg("cpd00239")))
    dat <- data.table(plotSpecActivity(sim, subs=paste0("EX_",dat.sub,"(e)"),useNames = F, rm_unused = T, ret_data = T))
    dat[,id:=str_extract(sub, "(?<=EX_).*?(?=\\(e\\))")]
    if(agora){
      spec.dic <- sapply(levels(dat$spec), function(s){trimws(agrep(gsub("_"," ",s), genome.desc$name, max=0.2, value=T))})
      names(spec.dic) <- levels(dat$spec)
      spec.dic$Salmonella_enterica_enterica_sv_Typhimurium_LT2 <- "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2"
      spec.dic$Desulfovibrio_desulfuricans_subsp_desulfuricans_DSM_642 <- "Desulfovibrio desulfuricans ND132"
      spec.dic$Enterococcus_faecalis_OG1RF_ATCC_47077 <- "Enterococcus faecalis EnGen0107 strain B594 chromosome"
      spec.dic$Lactobacillus_acidophilus_ATCC_4796 <- "Lactobacillus acidophilus strain DSM 20079 chromosome"
      spec.dic$Megasphaera_elsdenii_DSM_20460 <- "Megasphaera elsdenii 14-14"
      spec.dic$Veillonella_dispar_ATCC_17748 <- "Veillonella dispar strain NCTC11831 genome assembly"
      spec.dic$barkeri_iAF692 <- "barkeri_iAF692"
      dat$spec.name <- as.character(spec.dic[match(dat$spec, names(spec.dic))])
    }else{
      dat[,spec.name:=trimws(genome.desc$name[match(str_extract(spec, "GCF_[0-9]+\\.[0-9]"), gsub(".fna.gz","",genome.desc$file))])]  
    }
  }
  dat[,stat:=ifelse(mflux>0, "production", ifelse(mflux<0, "uptake", NA))]
  dat[spec=="barkeri_iAF692", spec.name:="Methanosarcina barkeri"]
}

# plot growth curve
growth.dt <- data.table()
growth.dt <- rbind(growth.dt, data.table(BacArena::plotGrowthCurve(sim.gapseq, use_biomass = F, ret_data = T))[, method:="gapseq"])
growth.dt <- rbind(growth.dt, data.table(BacArena::plotGrowthCurve(sim.modelseed, use_biomass = F, ret_data = T))[, method:="modelseed"])
growth.dt <- rbind(growth.dt, data.table(BacArena::plotGrowthCurve(sim.carveme, use_biomass = F, ret_data = T))[, method:="carveme"])
ggplot(growth.dt, aes(x=time, y=value)) + stat_summary(fun.y = "sum", geom="line") + facet_wrap(~method) + ylab("Number of organisms")


ex.gapseq    <- sim.activity(sim.gapseq, namespace.seed = T, sub.ferm=sub.ferm.seed)
ex.modelseed <- sim.activity(sim.modelseed, namespace.seed = T, sub.ferm=sub.ferm.seed)
ex.carveme   <- sim.activity(sim.carveme, namespace.seed = F, sub.ferm=sub.ferm.bigg)

dat.exp <- data.table(readODS::read_ods("dat/anaerobic_foodweb-reference.ods", na = c("NA", "<NA>"), range = "A1:U22", col_types = NA))
org.dic <- data.table(name=dat.exp$Organism)
org.dic[, name2:=sapply(name, function(org){ex.gapseq[spec.name %like% org, unique(spec.name)]})]
sub.dic <- data.table(name=c("Acetate", "Butyrate", "Ethanol", "Formate", "H2", "H2S", "Lactate", "Lactate", "Methane", "Propionate", "Succinate"),
                      seed=c("cpd00029", "cpd00211", "cpd00363", "cpd00047", "cpd11640", "cpd00239", "cpd00159", "cpd00221", "cpd01024", "cpd00141", "cpd00036"),
                      bigg=c("ac", "but", "etoh", "for", "h2", "h2s", "lac__L", "lac__D", "ch4", "ppa", "succ"),
                      vmh=c("ac", "but", "etoh", "for", "h2", "h2s", "lac_L", "lac_D", "ch4", "ppa", "succ"))

# get data from simulation files
extract_data <- function(t_gapseq, t_carveme, t_modelseed){
  dat.val <- data.table()
  for(i in 1:nrow(org.dic)){
    name1 <- org.dic$name[i]
    name2 <- unlist(org.dic$name2[i])
    if( length(name2)==0 ) next
    
    for(j in 2:ncol(dat.exp)){
      test.raw <- unlist(str_split(colnames(dat.exp)[j], "\\."))
      sub    <- test.raw[[1]]
      dir    <- test.raw[[2]]
      ref    <- dat.exp[Organism==name1][[j]]
      used   <- !ref==""
      sub.bigg <- sub.dic[name==sub, bigg]
      sub.seed <- sub.dic[name==sub, seed]
      sub.vmh  <- sub.dic[name==sub, vmh]
      flux.gapseq     <- ex.gapseq[spec.name==name2 & time==t_gapseq & id %in% sub.seed, mflux]; if(length(flux.gapseq)==0) flux.gapseq <- 0
      flux.modelseed  <- ex.modelseed[spec.name==name2 & time==t_modelseed & id %in% sub.seed, mflux]; if(length(flux.modelseed)==0) flux.modelseed <- 0
      flux.carveme    <- ex.carveme[spec.name==name2 & time==t_carveme & id %in% sub.bigg, mflux]; if(length(flux.carveme)==0) flux.carveme <- 0
      gapseq    <- ifelse( (dir=="prod" & flux.gapseq>0) | (dir=="up" & flux.gapseq<0), T, F)
      modelseed <- ifelse( (dir=="prod" & flux.modelseed>0) | (dir=="up" & flux.modelseed<0), T, F)
      carveme   <- ifelse( (dir=="prod" & flux.carveme>0) | (dir=="up" & flux.carveme<0), T, F)
      dat.val <- rbind(dat.val, data.table(org=name1, sub, dir, ref, used, gapseq, modelseed, carveme))
    }
  }
  
  #merge lactate_l and lactate_d
  dat.val.lac <- dat.val[sub=="Lactate", list(gapseq=as.logical(max(gapseq)), modelseed=as.logical(max(modelseed)), carveme=as.logical(max(carveme))), by=.(org,sub,dir,ref, used)]
  dat.val.unique <- rbind(dat.val[sub!="Lactate"], dat.val.lac)
  
  return(dat.val.unique)
}
dat.val <- extract_data(t_gapseq = 4, t_carveme = 4, t_modelseed = 6)

dat.val.melt <- melt(dat.val, id.vars = c("org", "sub", "dir", "ref", "used"), variable.name = "method", value.name = "prediction")
dat.val.melt[, val:=ifelse(used==T & prediction==T, "TP", ifelse(used==F & prediction==F, "TN", ifelse(used==T & prediction==F, "FN", "FP")))]
#dat.val.melt <- dat.val.melt[!(ref=="[6]" & val=="FN")] # remove FN from potential product list in Oliphant2019 (affect gapseq: 34, modelseed:50, carveme:51)

# validation
table(dat.val.melt[method=="gapseq", val]); round(.Last.value/dat.val.melt[method=="gapseq", .N],2)
table(dat.val.melt[method=="modelseed", val]); round(.Last.value/dat.val.melt[method=="modelseed", .N],2)
table(dat.val.melt[method=="carveme", val]); round(.Last.value/dat.val.melt[method=="carveme", .N],2)

