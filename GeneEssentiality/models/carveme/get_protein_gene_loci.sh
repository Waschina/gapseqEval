#!/bin/bash
for line in $(cat bsub_genes.lst) ; do
  printf "$line" >> bsub_genes_pos.tsv
  pos=`esearch -db gene -query "$line" | efetch -format docsum | xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop`
  printf "\t" >> bsub_genes_pos.tsv
  printf "$pos" >> bsub_genes_pos.tsv
  printf "\n" >> bsub_genes_pos.tsv
done



for line in $(cat ecol_genes.lst) ; do
  printf "$line" >> ecol_genes_pos.tsv
  pos=`esearch -db gene -query "$line" | efetch -format docsum | xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop`
  printf "\t" >> ecol_genes_pos.tsv
  printf "$pos" >> ecol_genes_pos.tsv
  printf "\n" >> ecol_genes_pos.tsv
done



for line in $(cat paer_genes.lst) ; do
  printf "$line" >> paer_genes_pos.tsv
  pos=`esearch -db gene -query "$line" | efetch -format docsum | xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop`
  printf "\t" >> paer_genes_pos.tsv
  printf "$pos" >> paer_genes_pos.tsv
  printf "\n" >> paer_genes_pos.tsv
done



for line in $(cat mgen_genes.lst) ; do
  printf "$line" >> mgen_genes_pos.tsv
  pos=`esearch -db gene -query "$line" | efetch -format docsum | xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop`
  printf "\t" >> mgen_genes_pos.tsv
  printf "$pos" >> mgen_genes_pos.tsv
  printf "\n" >> mgen_genes_pos.tsv
done



for line in $(cat sone_genes.lst) ; do
  printf "$line" >> sone_genes_pos.tsv
  pos=`esearch -db gene -query "$line" | efetch -format docsum | xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop`
  printf "\t" >> sone_genes_pos.tsv
  printf "$pos" >> sone_genes_pos.tsv
  printf "\n" >> sone_genes_pos.tsv
done
