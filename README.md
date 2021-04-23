# Structural basis for solubility in protein expression systems

![Twitter Follow](https://img.shields.io/twitter/follow/proteinqure?color=%23fffff&style=social)
![GitHub repo size](https://img.shields.io/github/repo-size/ProteinQure/cbh21-protein-solubility-challenge?style=plastic)

Large-scale protein production for biotechnology and biopharmaceutical applications rely on high protein solubility in expression systems. Solubility has been measured for a significant fraction of E. coli and S. cerevisiae proteomes and these datasets are routinely used to train predictors of protein solubility in different organisms. Thanks to continued advances in experimental structure-determination and modelling, many of these solubility measurements can now be paired with accurate structural models.

The challenge is mentored by [Christopher Ing](https://github.com/cing) and [Mark Fingerhuth](https://github.com/markf94).

## Challenge Aim

It is the objective of this project to use our provided dataset of protein structure and solubility value pairs in order to produce a solubility predictor with comparable accuracy to sequence-based predictors reported in the literature. The provided dataset to be used in this project is created by following the dataset curation procedure described in [the SOLart paper](https://academic.oup.com/bioinformatics/article/36/5/1445/5585748), and this hackathon project has a similar aim to this manuscript.

## The dataset

[more info here]

### Example Output
Your code should output a file called `predictions.csv` in the following format:

```
protein,solubility
P69829,83
P31133,62
```

whereby the `protein` column contains the Uniprot ID (corresponds to the filename of the PDB files) and the `solubility` column contains the predicted solubility value (can be `int` or `float`).

Note, that there are three (!) test subsets but you are expected to submit all the predictions in one file (not three) for the benchmarking system to work.

## Benchmarking System
The continuous integration script in `.github/workflows/ci.yml` will automatically build the `Dockerfile` on every commit to the `main` branch. This docker image will be published as your hackathon submission to `https://biolib.com/<YourTeam>/<TeamName>`. For this to work, make sure you set the `BIOLIB_TOKEN` and `BIOLIB_PROJECT_URI` accordingly as repository secrets. 

To read more about the benchmarking system [click here](https://www.notion.so/Benchmarking-System-46bfaeea0119490cb611688b493c589a).

## Say thanks

Give this repo a star ![GitHub Repo stars](https://img.shields.io/github/stars/ProteinQure/cbh21-protein-solubility-challenge?style=social)

Star the [ProteinQure](https://github.com/proteinqure) org on Github: ![GitHub Org's stars](https://img.shields.io/github/stars/ProteinQure?style=social)
