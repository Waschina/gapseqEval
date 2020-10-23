# Fermentation tests

All models can also be downloaded from ftp://ftp.rz.uni-kiel.de/pub/medsystbio/models/EnzymaticDataTestModels.zip

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
media_path="FT.csv" # Path might need to be adjusted
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

The script `extract_single_medium.R media/media.tsv [media_ID]` also exports a media file, that can directly be uploaded to the modelseed online platform. This media can than be used for modelseed's gapfilling.

Reconstruction itself is then performed on the webserver using the [ModelSEED-UI](https://modelseed.org/)