#!/bin/bash

ftp_paths=$(cat dat/protraits_genomes_quality.csv | tail -n +2 | cut -d '	' -f 4)

N=$(echo "$ftp_paths" | wc -l)
i=1
for path in $ftp_paths
do
    id=$(echo $path | grep -oE "GCA_[0-9]+")
    file1=$(ls *.fna.gz | grep $id)
    file2=$(ls *.faa.gz | grep $id)

    if [[ -n "$file1" ]]; then 
        #echo "$i/$N file exists $id (fna)"
        :
    else
        echo "$i/$N downloading $id (fna)"
        wget -nv `echo $path | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'` 
    fi 

    if [[ -n "$file2" ]]; then
        #echo "$i/$N file exists $id (faa)"
        :
    else
        echo "$i/$N downloading $id (faa)"
        wget -nv `echo $path | awk -F"/" '{print $0"/"$NF"_protein.faa.gz"}'` 
    fi

    i=$((i+1))
done
