library(data.table)

args <- commandArgs(trailingOnly=TRUE)
# $args[1] ... media bif table file
# $args[2] ... specific medium ID (e.g. "FT")

dt <- fread(args[1])
dt <- dt[medium == args[2]]

dt <- dt[,.(compounds = modelseed, name, maxFlux)]

dt.ms <- data.table(id = dt$compounds,
                    name = dt$name,
                    concentration = dt$maxFlux,
                    minflux = -1000,
                    maxflux = dt$maxFlux)

fwrite(dt, file = paste0(args[2],".csv"))
fwrite(dt.ms, file = paste0("media/",args[2],"_for_modelseed.csv"), quote = F, sep = "\t")