# Structural basis for solubility in protein expression systems

![GitHub repo size](https://img.shields.io/github/repo-size/ProteinQure/cbh21-protein-solubility-challenge?style=plastic)
![Twitter Follow](https://img.shields.io/twitter/follow/proteinqure?color=%23fffff&style=social)

Large-scale protein production for biotechnology and biopharmaceutical applications rely on high protein solubility in expression systems. Solubility has been measured for a significant fraction of E. coli and S. cerevisiae proteomes and these datasets are routinely used to train predictors of protein solubility in different organisms. Thanks to continued advances in experimental structure-determination and modelling, many of these solubility measurements can now be paired with accurate structural models.

The challenge is mentored by [Christopher Ing](https://github.com/cing) and [Mark Fingerhuth](https://github.com/markf94).

## Challenge Aim

It is the objective of this project to use our provided dataset of protein structure and solubility value pairs in order to produce a solubility predictor with comparable accuracy to sequence-based predictors reported in the literature. The provided dataset to be used in this project is created by following the dataset curation procedure of SOLart, and this hackathon project has a similar aim to this manuscript. 

### Example Output
You code should output a file called `predictions.csv` in the following format:

```
name,prediction
prot1,2
prot2,6
```

## Benchmarking System
The continuous integration script in `.github/workflows/ci.yml` will automatically build the `Dockerfile` on every commit to the `main` branch. This docker image will be published as your hackathon submission to `https://biolib.com/<YourTeam>/<TeamName>`. For this to work, make sure you set the `BIOLIB_TOKEN` and `BIOLIB_PROJECT_URI` accordingly as repository secrets. 

To read more about the benchmarking system [click here](https://www.notion.so/Benchmarking-System-46bfaeea0119490cb611688b493c589a).
