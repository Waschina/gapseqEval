# Fermentation tests

#### gapseq

Reconstructions with gapseq are calculated on the RZ-Cluster and CPLEX version 12.9.

The basic pipeline is:

```shell
model="MODELNAME" # same as genome file name without ".fna.gz" suffix

# Reaction & pathway prediction
$gapseq find -p all -b 200 -t [Bacteria / Archaea] $model.fna.gz

# Transporter prediction
$gapseq find-transport -b 200 $model.fna.gz

# Draft network reconstruction 
$gapseq draft -r $model-all-Reactions.tbl -t $model-Transporter.tbl -p $model-all-Pathways.tbl -u 200 -l 100 -c $model.fna.gz

# Gapfilling using FT-medium
media_path="/home/silvio/workspace/2018/gapseq/dat/media/FT.csv"
$gapseq fill -m $model-draft.RDS -n $media_path -c $model-rxnWeights.RDS -g $model-rxnXgenes.RDS -b 100 # for Archaea add: "-b archaea"
mv $model.RDS ${model}_FT.RDS
mv $model.xml ${model}_FT.xml

# Gapfilling using FT-medium w/o sugars
media_path="FT_nosugar.csv"
$gapseq fill -m $model-draft.RDS -n $media_path -c $model-rxnWeights.RDS -g $model-rxnXgenes.RDS -b 100
```



#### carveme

See script `carveme_recon.sh` for pipeline.



#### modelSEED

TODO