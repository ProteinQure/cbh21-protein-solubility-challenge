# Structural basis for solubility in protein expression systems

[![Twitter Follow](https://img.shields.io/twitter/follow/proteinqure?color=%23fffff&style=social)](https://twitter.com/proteinqure)
![GitHub repo size](https://img.shields.io/github/repo-size/ProteinQure/cbh21-protein-solubility-challenge?style=plastic)

Large-scale protein production for biotechnology and biopharmaceutical applications rely on high protein solubility in expression systems. Solubility has been measured for a significant fraction of E. coli and S. cerevisiae proteomes and these datasets are routinely used to train predictors of protein solubility in different organisms. Thanks to continued advances in experimental structure-determination and modelling, many of these solubility measurements can now be paired with accurate structural models.

The challenge is mentored by [Christopher Ing](https://github.com/cing) and [Mark Fingerhuth](https://github.com/markf94).

## Aim of the challenge

It is the objective of this project to use our provided dataset of protein structure and solubility value pairs in order to produce a solubility predictor with comparable accuracy to sequence-based predictors reported in the literature. The provided dataset to be used in this project is created by following the dataset curation procedure described in [the SOLart paper](https://academic.oup.com/bioinformatics/article/36/5/1445/5585748), and this hackathon project has a similar aim to this manuscript.

## The dataset

The process of generating the dataset is described in the SOLArt manuscript. At a high level, all experimentally tested E. coli and S. cerevisiae proteins were matched through Uniprot IDs to known crystallographic structures or high sequence similarity homology models. After balancing the fold types using CATH, a dataset containing a balanced spread of solubility values was produced. The resulting proteins for the training and testing of these models were prepared and disclosed in the supplemental material of this paper as a list of (Uniprot,PDB,Chain,Solubility) pairs. The PDB files were not included in this work so we had to re-extract them from [SWISS-MODEL](ttps://swissmodel.expasy.org/https://swissmodel.expasy.org/). Whenever a crystallographic structure was present, it was used, assuming high coverage over the Uniprot sequence. In some cases, the original PDB templates used within the original SOLArt paper had been superceded by improved templates, and we opted to take the highest resolution, highest sequence identity, models in our updated dataset. We stripped away all irrelevant chains and heteroatoms.

If issues are identified with individual structures, please refer to the Uniprot ID and manually investigate the best template. In some cases, we needed to improve structure correctness by modelling missing atoms/residues inside the Chemical Computing Group software MOE on a case-by-case basis.

The dataset can be found in the `data/` subdirectory - it is already divided into `training/` and `test/` data. The `training/` data comes with `solubility_values.csv` and `solublity_values.yaml` (same content just different format) which both contain the solubility target values for all the PDB files provided in that directory. Note that each PDB file is named after the Uniprot identifier of the respective protein and the `protein` column in the `solubility_values.csv` also contains the Uniprot identifiers.

The `test/` dataset consists of three different subdirectories (protein structures derived from different organisms and with different approaches) and you should **NOT** use them for any training. Only the `yeast_crystal_structs/` directory contains `solubility_values.csv` and `solublity_values.yaml` (same content just different format) files which you can use for some local testing & validation. In order to find out your performance on the entire test dataset you need to use the automated benchmarking system (see below).

### Example output
Your code should output a file called `predictions.csv` in the following format:

```
protein,solubility
P69829,83
P31133,62
```

whereby the `protein` column contains the Uniprot ID (corresponds to the filename of the PDB files) and the `solubility` column contains the predicted solubility value (can be `int` or `float`).

Note, that there are three (!) test subsets but you are expected to submit all the predictions in one file (not three) for the benchmarking system to work.

## Automated benchmarking system
The continuous integration script in `.github/workflows/ci.yml` will automatically build the `Dockerfile` on every commit to the `main` branch. This docker image will be published as your hackathon submission to `https://biolib.com/<YourTeam>/<TeamName>`. For this to work, make sure you set the `BIOLIB_TOKEN` and `BIOLIB_PROJECT_URI` accordingly as repository secrets. 

To read more about the benchmarking system [click here](https://www.notion.so/Benchmarking-System-46bfaeea0119490cb611688b493c589a).

## Say thanks

Give this repo a star: ![GitHub Repo stars](https://img.shields.io/github/stars/ProteinQure/cbh21-protein-solubility-challenge?style=social)

Star the [ProteinQure](https://github.com/proteinqure) org on Github: ![GitHub Org's stars](https://img.shields.io/github/stars/ProteinQure?style=social)
