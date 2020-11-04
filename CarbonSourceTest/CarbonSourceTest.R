library(data.table)
library(sybil)
library(ggplot2)
library(reshape2)
library(stringr)

cs.sub <- fread("dat/sub2pwy.csv") # gapseq file containing carbon sources + exchange ids

# file from: http://protraits.irb.hr/data.html
protraits.db <- fread("dat/ProTraits_binaryIntegratedPr0.95.txt")
cs.src <- colnames(protraits.db)[3:111]
protraits.cs <- protraits.db[protraits.db[, Reduce(`|`, lapply(.SD, `!=`, "?")),.SDcols = cs.src],] # organisms with at least one carbon source prediction

#
# load models (download from ftp://ftp.rz.uni-kiel.de/pub/medsystbio/models/CarbonSourceTestModels.zip)
#
modelseed <- readRDS("models/protraits_modelseed2.RDS")
carveme   <- readRDS("models/protraits_carveme2.RDS")
gapseq <- readRDS("models/protraits.old_20200910.RDS")
gapseq.id <- gsub("_genomic","",names(gapseq))
modelseed.id <- gsub("_genomic","",gsub(".fna.sbml$","",names(modelseed)))
carveme.id <- gsub("_genomic","",gsub("_protein.faa","", names(carveme)))

#
# matching of carbon source names
#
dt.cs  <- data.table(name=gsub("_"," ",cs.src))
idx    <- match(dt.cs$name, tolower(cs.sub$name))
idx[is.na(idx)] <- match( dt.cs$name[is.na(idx)], tolower(cs.sub$altname))
idx[is.na(idx)] <- match( dt.cs$name[is.na(idx)], tolower(str_remove(cs.sub$name, "^.-")))
idx[is.na(idx)] <- match( str_remove(dt.cs$name[is.na(idx)],"^.-"), tolower(cs.sub$name))
idx[is.na(idx)] <- match( str_remove(dt.cs$name[is.na(idx)],"^.-"), tolower(cs.sub$altname))
dt.cs[,seed.name:=cs.sub$name[idx]]
dt.cs[,seed.ex:=cs.sub$exid_seed[idx]]
dt.cs[,vmh.ex:=cs.sub$exid[idx]]
dt.cs <- dt.cs[!is.na(seed.ex) & !seed.ex=="" & !is.na(vmh.ex) & !vmh.ex==""]
# manual modifications
dt.cs[name=="dextrin", seed.ex:="EX_cpd11976_e0"] # dextrin not really distingushable from maltodextrin and dextrin=maltodextrin for seed db



# refseq genomes for protraits organism
tax2genome <- fread("dat/protraits_genomes.csv", fill = T, header=F, sep="\t"); tax2genome$V4 <- str_extract(tax2genome$V3, "GCA_.*")

set_diet <- function(model, medium, uptake=-100,verbose=TRUE){
  ex <- findExchReact(model)
  ex <- ex[grep("^EX_",ex@react_id),] # reduce to real exchange reactions only
  if(is.null(medium)){
    medium <- ex@react_id[which(ex@lowbnd < 0)]
    uptake <- ex@lowbnd[which(ex@lowbnd < 0)]
  } 
  
  ub <- model@uppbnd
  lb <- model@lowbnd
  if(length(uptake) == 1)
    uptake <- rep(uptake, length(medium))
  
  lb[ex@react_pos] <- 0
  idx <- match(medium, react_id(model))
  if(sum(is.na(idx))>0 & verbose) cat("Not found in model:", medium[is.na(idx)], "\n")
  lb[na.omit(idx)] <- uptake[!is.na(idx)]
  model.constrainted <- changeBounds(model, react=1:length(lb), ub=ub, lb=lb)
  
  return(model.constrainted)
}

csource2 <- function(mod, cs, medium=NULL, verbose=FALSE, esource=F, GrowthExNa=0, db="seed"){
  if( sybil::SYBIL_SETTINGS("SOLVER") == "cplexAPI" ) solver_ok=1 else if( sybil::SYBIL_SETTINGS("SOLVER") == "glpkAPI" ) solver_ok=5
  
  if( length(medium) == 0 ){
    if ( verbose ) print("No medium specified, using glucose minimal medium!")
    minmedia <- fread("dat/MM_glu.csv")
    medium <- paste0("EX_",minmedia[,compounds],"_e0")
  }
  
  if( esource ){
    if( db=="seed"){
      mql <- "cpd15499[c0]"; mqn  <- "cpd15500[c0]"
      uql <- "cpd15561[c0]"; uqn  <- "cpd15560[c0]"
      nadh<- "cpd00004[c0]"; nad  <- "cpd00003[c0]"
      h   <- "cpd00067[c0]"
    }else{
      mql <- "mql8[c]"; mqn   <- "mqn8[c]"
      uql <- "u8h2[c]"; uqn   <- "u8[c]"
      nadh<- "nadh[c]"; nad   <- "nad[c]"
      h   <- "h[c]"
    }
    mod <- addReact(mod, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if ubiquinone pool can be recycled 
    mod <- addReact(mod, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if menaquinone pool can be recycled 
    mod <- addReact(mod, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1),lb=0, ub=1000) # check if nadh can be recycled
    mod <- changeObjFunc(mod, react=c("ESP1", "ESP2", "ESP3"), obj_coef=c(1,1,1))
  }
  med_withoutCS <- setdiff(medium, "EX_cpd00027_e0") # remove glucose
  if( db=="vmh"){
    dic <- fread("dat/SEED2VMH_translation_edited.csv", header=F)
    idx <- match(gsub("_e0","\\(e\\)",med_withoutCS), dic$V1)
    med_withoutCS <- dic$V2[idx]
  }
  dfcs <- data.frame()
  for( carbon in cs){
    if( !carbon %in% react_id(mod) ){
      dfcs <- rbind(dfcs, data.frame(csource=carbon, usage=FALSE, growth=GrowthExNa)) # if no exchange is avalable set set to default value
      next
    }
    med <- c(med_withoutCS, carbon)
    model <- set_diet(mod, med, verbose=verbose)
    sol <- sybil::optimizeProb(model, retOptSol=F)
    usage = sol$stat==solver_ok & round(sol$obj,3)>0
    growth= ifelse(sol$stat==solver_ok, sol$obj, 0)
    dfcs <- rbind(dfcs, data.frame(csource=carbon, usage, growth))
    
    if( verbose & db == "seed" ){
      idx.active <- which(sol$fluxes != 0)
      mod@react_name[idx.active]  
      idx.active.met <- idx.active[which(colSums(abs(mod@S[grep(str_extract(carbon, "cpd[0-9]+"), mod@met_id), idx.active]))!=0)]
      mod@react_name[idx.active.met]
      mod@react_id[idx.active.met]      
    }
  }
  return(dfcs)
}

library(foreach)
library(doParallel)
library(cplexAPI)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
registerDoParallel(3) # 3 cores; detectCores()
dt.cs.predict <- foreach(i=1:nrow(protraits.cs), .combine=rbind) %dopar%{
  cat("\r",i,"/",nrow(protraits.cs))
  tax.id <- protraits.cs$Tax_ID[i]
  genome.id <- tax2genome[V1==tax.id,V4]
  if( !is.na(genome.id) ){
    idx.gapseq     <- match(genome.id, gsub("GCF","GCA", gapseq.id)) # for old genome files
    idx.modelseed <- match(genome.id, gsub("GCF","GCA", modelseed.id))
    idx.carveme   <- match(genome.id, gsub("GCF","GCA", carveme.id)) 
    if( !all(is.na(idx.gapseq)) & !all(is.na(idx.modelseed)) & !all(is.na(idx.carveme)) ){
      mod.gapseq    <- gapseq[[   na.omit(idx.gapseq)[1] ]]
      mod.modelseed <- modelseed[[na.omit(idx.modelseed)[1] ]]
      mod.carveme   <- carveme[[  na.omit(idx.carveme)[1] ]]
      
      cs.pot <- protraits.cs[Tax_ID==tax.id, c(1:111)]
      cs <- data.table::melt(cs.pot, id.vars = c("Organism_name", "Tax_ID"))[value!="?" & !is.na(value)]
      cs.sel <- dt.cs[name %in% cs$variable]  
      if( nrow(cs.sel)>0 ){
        growth.gapseq    <- csource2(mod.gapseq   , cs=cs.sel$seed.ex,esource=TRUE)
        growth.modelseed <- csource2(mod.modelseed, cs=cs.sel$seed.ex,esource=TRUE)
        growth.carveme   <- csource2(mod.carveme  , cs=cs.sel$vmh.ex,esource=TRUE, db="vmh")
        
        #cs.sel[,org:= genome.id]; 
        cs.sel[,org:= tax.id]; cs.sel[,org.name:= protraits.cs$Organism_name[i]]
        cs.sel[,org.gapseq:= gapseq.id[na.omit(idx.gapseq)[1] ]]
        cs.sel[,protraits:=cs[match(cs.sel$name, cs$variable), value==1]]
        cs.sel[,gapseq:=growth.gapseq$usage]
        cs.sel[,modelseed:=growth.modelseed$usage]
        cs.sel[,carveme:=growth.carveme$usage]
        cs.sel
      }
    } else data.table()
  } else data.table()
}

#dt.cs.predict <- readRDS("dat/protraits.old_20200910.RDS")

dt.cs.predict <- dt.cs.predict[name != "d-lyxose"] # 500 TN for all methods (only negative usage in database, would introduce big bias in figure)

dt.cs.predict.melt <- data.table::melt(dt.cs.predict[,.(org,protraits,gapseq,modelseed,carveme,name,seed.ex)], id.vars = c("org","protraits", "name", "seed.ex"), variable.name = "method")
dt.cs.predict.melt[value==T & protraits==T, validation:="TP"]
dt.cs.predict.melt[value==F & protraits==F, validation:="TN"]
dt.cs.predict.melt[value==T & protraits==F, validation:="FP"]
dt.cs.predict.melt[value==F & protraits==T, validation:="FN"]

table(dt.cs.predict.melt[method=="carveme",validation]); round(.Last.value / dt.cs.predict.melt[method=="carveme",.N], 2) # 0.36 0.11 0.29 0.24
table(dt.cs.predict.melt[method=="modelseed",validation]); round(.Last.value / dt.cs.predict.melt[method=="modelseed",.N], 2) # 0.29 0.05 0.35 0.31
table(dt.cs.predict.melt[method=="gapseq",validation]); round(.Last.value / dt.cs.predict.melt[method=="gapseq",.N], 2) # 0.13 0.11 0.29 0.47

dt.cs.predict.plot <- dt.cs.predict.melt[,.N, by=.(method, validation)]
dt.cs.predict.plot$method <- factor(dt.cs.predict.plot$method, levels = c("carveme", "gapseq", "modelseed"), labels= c("CarveMe", "gapseq", "ModelSEED"), ordered = T)
dt.cs.predict.plot$validation <- factor(dt.cs.predict.plot$validation, levels = c("FN", "FP", "TN", "TP"),labels = c("False negative", "False positive", "True negative", "True positive"), ordered = T)
dt.cs.predict.plot$validation <- factor(dt.cs.predict.plot$validation, levels = c("True positive", "True negative", "False positive", "False negative"))
ggplot(data = dt.cs.predict.plot) + geom_bar(stat="identity",aes(x=validation, y=N,fill=method),position="dodge2", color="black") + 
  #geom_text(aes(x=validation, y=N, label=N), position = position_dodge2(width = 1)) + # wrong order?????
  xlab("") + ylab("Carbon source comparison") + labs(fill="") + theme_bw(base_size = 14) + 
  scale_y_continuous(limits=c(0, 875), expand = c(0, 0)) + #scale_x_continuous(expand = c(0, 0)) + 
  #scale_fill_discrete() +
  scale_fill_manual(values=c("#377eb8", "#e41a1c", "#FFB200")) +
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank())
#ggsave("~/uni/gapseq_doc/img/protraits.pdf", width=6, height=5)


# good/bad predicted carbon sources
dt.quali.all <- dt.cs.predict.melt[,list(FN=sum(validation=="FN"), FP=sum(validation=="FP"),
                                         TN=sum(validation=="TN"), TP=sum(validation=="TP")), by=name]
dt.quali.all[,accuracy:=(TN+TP)/(FN+FP+TN+TP)]
dt.quali.all[,sensitivity:=(TP)/(TP+FN)]
dt.quali.all[,specificity:=(TN)/(TN+FP)]
dt.quali.all[,N:=TN+TP+FP+FN]
dt.quali.all[order(accuracy)][FN+FP+TN+TP>100]
