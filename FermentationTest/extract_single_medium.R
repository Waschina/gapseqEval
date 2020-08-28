library(data.table)

args <- commandArgs(trailingOnly=TRUE)
# $args[1] ... media bif table file
# $args[2] ... specific medium ID (e.g. "FT")

dt <- fread(args[1])
dt <- dt[medium == args[2]]

dt <- dt[,.(compunds = modelseed, name, maxFlux)]

fwrite(dt, file = paste0(args[2],".csv"))