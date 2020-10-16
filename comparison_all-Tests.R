library(readODS)
library(data.table)
library(ggplot2)
library(reshape2)

getwd()
dat.file <- "doc/dat/summary_tests.ods"
dat.file <- "summary_tests.ods" # please adjust

confusion.eval <- function(df){
  fp <- as.numeric(df[1,2])
  fn <- as.numeric(df[2,1])
  tp <- as.numeric(df[2,2])
  tn <- as.numeric(df[1,1])
  return(data.table(sensitivity=tp/(tp+fn), specificity=tn/(tn+fp), similarity=(tp+tn)/(tp+tn+fp+fn)))
}

read_dat <- function(file, range, test.name){
  df <- read_ods(file, range=range, col_names = F)
  df <- df[,colSums(is.na(df))<nrow(df)] # remove NA columns
  df.carveme   <- cbind(confusion.eval(df[1:2,2:3]), data.table(Method="CarveMe"))
  df.gapseq    <- cbind(confusion.eval(df[1:2,4:5]), data.table(Method="gapseq"))
  df.modelseed <- cbind(confusion.eval(df[1:2,6:7]), data.table(Method="ModelSEED"))
  
  dat.dt <- rbind(df.carveme, df.gapseq, df.modelseed)
  dat.dt$test <- test.name
  return(dat.dt)
}

test.ferm <- read_dat(dat.file, range="C4:K5", test="Fermentation products")
test.ess  <- read_dat(dat.file, range="C13:K14", test="Gene essentiality")
test.enz  <- read_dat(dat.file, range="C22:K23", test="Enzyme activity")
test.esrc <- read_dat(dat.file, range="C31:K32", test="Energy source")
test.fdwb <- read_dat(dat.file, range="C40:K41", test="Anaerobic food web")

test.dt <- rbind(test.ferm, test.ess, test.enz, test.esrc, test.fdwb)

test.dt[,list(similarity=mean(similarity), sensitivity=mean(sensitivity), specificity=mean(specificity)),by=Method]

#ggplot(test.dt, aes(x=specificity, y=sensitivity, shape=test, fill=Method)) +
#  scale_fill_manual(values=c("#377eb8", "#e41a1c", "#FFB200")) +
#  scale_shape_manual(values = 21:24) +
  #geom_abline(slope = 1, intercept = 0,
  #            na.rm = FALSE, show.legend = NA) +
  #geom_point(size = 2.5) +
#  geom_point(size = 4) +
#  theme_bw() + theme(legend.title=element_blank(), 
#                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#  scale_x_continuous(limits = c(0.55,1)) + scale_y_continuous(limits = c(0,1)) +
#  guides(fill = guide_legend(override.aes=list(shape=NA)))
#ggsave("doc/img/comparison_all-Tests.pdf", width=4.75, height=3)
#ggsave("~/Documents/comparison_all-Tests.pdf", width=4.75, height=3)

ggplot(test.dt, aes(x=specificity, y=sensitivity, fill=Method)) +
  scale_fill_manual(values=c("#377eb8", "#e41a1c", "#FFB200")) +
  geom_point(size = 4, shape = 21) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0.55,1)) + scale_y_continuous(limits = c(0,1)) +
  facet_wrap("test")
ggsave("~/Documents/comparison_all-Tests.pdf", width=5, height=3.8)




