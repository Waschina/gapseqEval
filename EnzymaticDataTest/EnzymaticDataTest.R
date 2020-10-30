library(data.table)
library(stringr)
library(foreach)
library(doParallel)
library(sybil)
library(ggplot2)

ec_seed     <- fread("dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv")
ec_vmh      <- fread("dat/vmh_reactions.csv")
org.dt <- fread("dat/bacdive_genomes.log")
db.dat  <- readRDS("dat/type-strains_data.RDS")


#
# 1) Enzymatic data test
#
# carveme, modelseed, and gapseq models (3.5gb) have to be downloaded from: 
# ftp://ftp.rz.uni-kiel.de/pub/medsystbio/models/EnzymaticDataTestModels.zip
#
# Loading the model files is a bit memory intensive ~50gb
# Alternatively, the result file is also provided in this repository and 
# evaluation can be continued in step 2), see below

# download EnzymaticDataTestModels.zip and extract files to models
models.gapseq    <- readRDS("models/bacdive-20200929.RDS")
models.carveme   <- readRDS("models/bacdive-carveme_20200130.RDS")
models.modelseed <- readRDS("models/bacdive-modelseed_20200213.RDS")

orgmatch.dt <- org.dt[!is.na(V3)]

models.gapseq.id <- str_extract(names(models.gapseq), "GCF_[0-9]+.[0-9]")
models.modelseed.id <- str_extract(names(models.modelseed), "GCF_[0-9]+.[0-9]")
models.carveme.id <- str_extract(names(models.carveme), "GCF_[0-9]+.[0-9]")

registerDoParallel(1) # potentially high memory consumption, doesn't take long with one core though
val.dt <- foreach(i=1:nrow(orgmatch.dt), .combine=rbind) %dopar%{
  #i <- 33
  bacdive_nr <- orgmatch.dt[i,V6]
  bacdive_exp<- db.dat[[match(bacdive_nr, names(db.dat))]]
  if( length(bacdive_exp$morphology_physiology) > 0 ){
    if( "enzymes" %in% names(bacdive_exp$morphology_physiology) ){
      if( orgmatch.dt[i,V2] %in% models.gapseq.id && orgmatch.dt[i,V2] %in% models.modelseed.id && orgmatch.dt[i,V2] %in% models.carveme.id ){
        bacdive_enzyme <- bacdive_exp$morphology_physiology$enzymes
        
        mod.gapseq    <- models.gapseq[[ match(orgmatch.dt[i,V2], models.gapseq.id) ]]
        mod.modelseed <- models.modelseed[[ match(orgmatch.dt[i,V2], models.modelseed.id) ]]
        mod.carveme   <- models.carveme[[ match(orgmatch.dt[i,V2], models.carveme.id) ]]
        
        mod.gapseq.rxn <- str_extract(mod.gapseq@react_id, "rxn[0-9]+")
        mod.modelseed.rxn <- str_extract(mod.modelseed@react_id, "rxn[0-9]+")
        mod.carveme.rxn  <- mod.carveme@react_id
        
        pred.tmp <- data.table(org=orgmatch.dt[i,V1], id=bacdive_nr, genome=orgmatch.dt[i,V2], ec=bacdive_enzyme$ec_number, activity=bacdive_enzyme$activity, enzyme=bacdive_enzyme$enzyme)
        
        pred.tmp[,gapseq:=any(unlist(str_split(ec_seed[`External ID`==ec,`MS ID`],"\\|")) %in% mod.gapseq.rxn), by=ec]
        pred.tmp[,modelseed:=  any(unlist(str_split(ec_seed[`External ID`==ec,`MS ID`],"\\|")) %in% mod.modelseed.rxn), by=ec]
        pred.tmp[,carveme:= any(ec_vmh[ecnumber==ec,abbreviation] %in% mod.carveme.rxn), by=ec]
        
        pred.tmp        
        
      } else data.table()
    } else data.table()
  } else data.table()
}
stopImplicitCluster()


#val.dt <- readRDS("dat/bacdive_enzyme-comp_20200929.RDS") # run used in manuscript

val.dt <- val.dt[ec!="" & activity %in% c("+", "-")] # only unique data & with EC number
val.dt <- val.dt[!ec %in% c("3.1.3.1", "3.1.3.2")] # removed because they are catalyzing the same reactions (kegg, metacyc)

# consider only those EC numbers for which matching exists
mis.vmh <- sapply(unique(val.dt$ec), function(ec){ec_vmh[ecnumber==ec,abbreviation]})
mis.seed<- sapply(unique(val.dt$ec), function(ec){ec_seed[`External ID`==ec,`MS ID`]})
mis.ec_match <- unique(c(names(mis.vmh[which(sapply(mis.vmh,length)==0)]), names(mis.seed[which(sapply(mis.seed,length)==0)])))
val.dt <- val.dt[!ec %in% mis.ec_match]

#
# 2) Enzymatic data test evaluation
#

# figure
val.dt.melt <- melt(val.dt, id.vars = c("org","ec", "activity", "enzyme", "id", "genome"), variable.name = "method")
val.dt.melt[activity=="+" & value==TRUE, validation:="True positive"]
val.dt.melt[activity=="-" & value==FALSE, validation:="True negative"]
val.dt.melt[activity=="+" & value==FALSE, validation:="False negative"]                     
val.dt.melt[activity=="-" & value==TRUE, validation:="False positive"]                     
val.dt.table <- as.data.table(table(val.dt.melt[,.(method,validation)]))
val.dt.table$method <- factor(val.dt.table$method)
val.dt.table$validation <- factor(val.dt.table$validation, levels = c("True positive", "True negative", "False positive", "False negative"))
ggplot(data = val.dt.table) + geom_bar(stat="identity",aes(x=validation, y=N,fill=as.character(method)),position="dodge2", color="black") + 
  xlab("") + ylab("Enzyme activity comparison") + labs(fill="") + theme_bw(base_size = 14) + 
  scale_y_continuous(limits=c(0, 6000), expand = c(0, 0)) + #scale_x_continuous(expand = c(0, 0)) + 
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values=c("#377eb8", "#e41a1c", "#FFB200"))
#ggsave("enzyme-test.pdf", width=6, height=5)  



#
# 3) Independent and identical sampling (check for bias)
#

val.dt[,.N,by=ec][order(N)]
ec.high.cov <- val.dt[,.N,by=ec][N>100, ec] # consider EC numbers with more than 100 tests
B=500 # bootstrap sample size
resampling.dt <- data.table()
for(i in 1:B){
  val.dt.rand <- data.table()
  for(e in ec.high.cov){
    val.dt.rand <- rbind(val.dt.rand, val.dt[ec == e][sample(.N, 100)])
  }
  val.dt.rand.melt <- melt(val.dt.rand, id.vars = c("org","ec", "activity", "enzyme", "id", "genome"), variable.name = "method")
  val.dt.rand.melt[activity=="+" & value==TRUE, validation:="True positive"]
  val.dt.rand.melt[activity=="-" & value==FALSE, validation:="True negative"]
  val.dt.rand.melt[activity=="+" & value==FALSE, validation:="False negative"]                     
  val.dt.rand.melt[activity=="-" & value==TRUE, validation:="False positive"]                     
  val.dt.rand.table <- as.data.table(table(val.dt.rand.melt[,.(method,validation)]))
  
  gapseq.sensitivity <- val.dt.rand.table[method=="gapseq" & validation=="True positive", N] / (val.dt.rand.table[method=="gapseq" & validation=="True positive", N] + val.dt.rand.table[method=="gapseq" & validation=="False negative", N])
  gapseq.specificity <- val.dt.rand.table[method=="gapseq" & validation=="True negative", N] / (val.dt.rand.table[method=="gapseq" & validation=="True negative", N] + val.dt.rand.table[method=="gapseq" & validation=="False positive", N])
  carveme.sensitivity <- val.dt.rand.table[method=="carveme" & validation=="True positive", N] / (val.dt.rand.table[method=="carveme" & validation=="True positive", N] + val.dt.rand.table[method=="carveme" & validation=="False negative", N])
  carveme.specificity <- val.dt.rand.table[method=="carveme" & validation=="True negative", N] / (val.dt.rand.table[method=="carveme" & validation=="True negative", N] + val.dt.rand.table[method=="carveme" & validation=="False positive", N])
  modelseed.sensitivity <- val.dt.rand.table[method=="modelseed" & validation=="True positive", N] / (val.dt.rand.table[method=="modelseed" & validation=="True positive", N] + val.dt.rand.table[method=="modelseed" & validation=="False negative", N])
  modelseed.specificity <- val.dt.rand.table[method=="modelseed" & validation=="True negative", N] / (val.dt.rand.table[method=="modelseed" & validation=="True negative", N] + val.dt.rand.table[method=="modelseed" & validation=="False positive", N])
  
  resampling.dt <- rbind(resampling.dt, data.table(nr=i, gapseq.sensitivity, gapseq.specificity, carveme.sensitivity, carveme.specificity, modelseed.sensitivity, modelseed.specificity))
}

alldat <- rbindlist(list(data.table(method="gapseq",statistic="sensitivity", value=val.dt.table[method=="gapseq" & validation=="True positive", N] / (val.dt.table[method=="gapseq" & validation=="True positive", N] + val.dt.table[method=="gapseq" & validation=="False negative", N])),
                         data.table(method="gapseq",    statistic="specificity", value=val.dt.table[method=="gapseq" & validation=="True negative", N] / (val.dt.table[method=="gapseq" & validation=="True negative", N] + val.dt.table[method=="gapseq" & validation=="False positive", N])),
                         data.table(method="carveme",   statistic="sensitivity", value=val.dt.table[method=="carveme" & validation=="True positive", N] / (val.dt.table[method=="carveme" & validation=="True positive", N] + val.dt.table[method=="carveme" & validation=="False negative", N])),
                         data.table(method="carveme",   statistic="specificity", value=val.dt.table[method=="carveme" & validation=="True negative", N] / (val.dt.table[method=="carveme" & validation=="True negative", N] + val.dt.table[method=="carveme" & validation=="False positive", N])),
                         data.table(method="modelseed", statistic="sensitivity", value=val.dt.table[method=="modelseed" & validation=="True positive", N] / (val.dt.table[method=="modelseed" & validation=="True positive", N] + val.dt.table[method=="modelseed" & validation=="False negative", N])),
                         data.table(method="modelseed", statistic="specificity", value=val.dt.table[method=="modelseed" & validation=="True negative", N] / (val.dt.table[method=="modelseed" & validation=="True negative", N] + val.dt.table[method=="modelseed" & validation=="False positive", N]))))

resampling.melt <- rbindlist(list(data.table(method="gapseq", statistic="sensitivity",    dat=resampling.dt$gapseq.sensitivity),
                                  data.table(method="gapseq", statistic="specificity",    dat=resampling.dt$gapseq.specificity),
                                  data.table(method="carveme", statistic="sensitivity",   dat=resampling.dt$carveme.sensitivity),
                                  data.table(method="carveme", statistic="specificity",   dat=resampling.dt$carveme.specificity),
                                  data.table(method="modelseed", statistic="sensitivity", dat=resampling.dt$modelseed.sensitivity),
                                  data.table(method="modelseed", statistic="specificity", dat=resampling.dt$modelseed.specificity)
)
)
resampling.melt$method <- factor(resampling.melt$method)
ggplot(resampling.melt, aes(x=method, y=dat, fill=method)) + 
  geom_boxplot(color="black") +
  geom_point(data = alldat, aes(x=method, y=value, color=method), shape=18, size=4) +   
  xlab("") + ylab("") + labs(fill="Sampling", color="All data") +
  facet_wrap(~statistic) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) + # needed to overwrite shape
  theme_bw(base_size = 14) + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values=c("#377eb8", "#e41a1c", "#FFB200")) + scale_color_manual(values=c("#377eb8", "#e41a1c", "#FFB200"))
#ggsave("enzyme-test_sampling.pdf", width=5, height=4)


