calc_seq_overlap <- function(astart, aend, bstart, bend, a.contig = "", b.contig = rep("", length(bstart))) {
  out <- numeric(length = length(bstart))
  for(i in seq_along(bstart)){
    if(a.contig == b.contig[i])
      out[i] <- length(intersect(astart:aend,bstart[i]:bend[i]))/((length(astart:aend)+length(bstart[i]:bend[i]))/2)
    else
      out[i] <- 0
  }
  return(out)
}