# Example Challenge

This repository is a template for creating challenges. You can change this README to describe your challenge in detail.

The participants will fork your challenge repository at the beginning of the hackathon.

## Challenge Aim

The aim of this challenge is to predict adenine amino acid count. You are given a FASTA file using the `--infile` argument. For each protein, output your prediction for the number adenine amino acids in that protein. 

### Example Output
You code should output a file called `predictions.csv` in the following format:

```
name,prediction
prot1,2
prot2,6
```

## Benchmarking System
The continuous integration script in `.github/workflows/ci.yml` will automatically build the `Dockerfile` on every commit to the `main` branch. This docker image will be published as your hackathon submission to `https://biolib.com/<YourTeam>/<TeamName>`. For this to work, make sure you set the `BIOLIB_TOKEN` and `BIOLIB_PROJECT_URI` accordingly as repository secrets. 
