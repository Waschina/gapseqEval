prepare_MS_gene_table <- function(ms.gene.file, mod) {
  require(data.table)
  require(jsonlite)
  
  msg <- jsonlite::fromJSON(ms.gene.file)
  msg <- data.table(msg$features)
  msg <- msg[,.(id, location)]
  
  msg$contig  <- unlist(lapply(msg$location, FUN = function(x) return(x[1])))
  msg$sstart  <- unlist(lapply(msg$location, FUN = function(x) return(x[2])))
  msg$sdir    <- unlist(lapply(msg$location, FUN = function(x) return(x[3])))
  msg$slength <- unlist(lapply(msg$location, FUN = function(x) return(x[4])))
  
  msg$sstart  <- as.numeric(msg$sstart)
  msg$slength <- as.numeric(msg$slength)
  
  msg[sdir == "+", send := sstart + slength]
  msg[sdir == "-", send := sstart - slength]
  
  msg[, id := gsub("fig|", "", id, fixed = T)]
  
  msg <- msg[,.(mod.gene = id, contig, start = sstart, end = send)]
  
  msg <- msg[mod.gene %in% mod@allGenes]
  
  return(msg)
}
