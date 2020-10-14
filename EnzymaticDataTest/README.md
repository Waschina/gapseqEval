# Enzymatic data tests

- ``EnzymaticDataTest.R`` contains code to redo analysis and create figures that were used in the manuscript
- carveme, modelseed and gapseq models are available at ftp://ftp.rz.uni-kiel.de/pub/medsystbio/models/EnzymaticDataTestModels.zip


## genome download
The refseq ids for BacDive reference strains, which were downloaded and used for metabolic model construction, are listed in ``EnzymaticDataTest/dat/bacdive_genomes.log``.


## gapseq

Reconstructions with gapseq were calculated on the RZ-Cluster and CPLEX version 12.9.

The basic pipeline is:

```shell
base=$(basename "$file")
if [[ $base == *.gz ]]; then
    id_tmp=${base%.*}
    id=${id_tmp%.*}
else
    id=${base%.*}
fi

gapseq find -v 0 -k -b 200 -p all -m bacteria "$file"

gapseq find-transport -k -b 200 "$file"

gapseq draft -r "$gapseq_files/$pwy_out/$id-all-Reactions.tbl" -t "$gapseq_files/$pwy_out/$id-Transporter.tbl" -c "$file" -u 200 -l 100 -a 1 -p "$gapseq_files/$pwy_out/$id-all-Pathways.tbl"

gapseq fill -m "$gapseq_files/out/${id}-draft.RDS" -n ./gapseq/dat/media/ALLmed.csv -c "$gapseq_files/out.tmp/${id}-rxnWeights.RDS" -b 100 -g "$gapseq_files/out.tmp/${id}-rxnXgenes.RDS"
```

## carveme

parallel --eta carve -r {} -o ./ ::: *.faa.gz


## modelSEED

Reconstruction itself was performed on the webserver using the [ModelSEED-UI](https://modelseed.org/).
Selenium in python3 was used for automatisation.
