#!/bin/bash
cd genomes.faa

cm_recon() {
  
  gunzip -k ${1}_protein.faa.gz
  
  ~/.local/bin/./carve "${1}_protein.faa" --gapfill $4 --mediadb ../media/media.tsv -o "../models/carveme/$1.xml" --diamond-args="-b6"
  mv "${1}_protein.tsv" "../models/carveme/$1.tsv"
  
  rm ${1}_protein.faa
  
  echo "DONE(carveme): $1 - $4"
}

#export python
export -f cm_recon


cat ../organisms2.csv | \
  parallel --colsep '\t' --header : cm_recon $id $media

