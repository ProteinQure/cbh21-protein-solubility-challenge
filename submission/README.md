# Bits Please

This is our submission for the CBH2021 hackathon. 

Our team is from Copenhagen and we think proteins are cool 🙌.

## Project Description

Our team aimed to create a model to predict protein solubility. We used a random forest with a selected assortment of features. 

## Team Members

    Javier Marchena
    Roosa Varjus
    Caroline Linnea Elin Lennartsson
    Marida Ianni-Ravn
    Henrietta Holze

## Features 
We chose a selected amount of features from the article "SOLart: a structure-based method to predict protein
solubility and aggregation" by Qingzhen Hou1, Jean Marc Kwasigroch,  Marianne Rooman and
Fabrizio Pucci.  In addition, we also added some new features to investigate the properties of solubility further. 

The features we investigated in described below; 

* Accessible surface area in Å^2. 
* The length of the amino acid sequence. 
* Fraction of the surface divided by sequence length. 
* Fraction of moderatly buried beta residues. 
* Fraction of moderatly buried alfa residues. 
* Fraction of exposed buried alfa residues. 
* Fraction of K residues minus the fraction of R residues.   
* Fraction of negativly charged residues. 
* Fraction of positivly charged residues.  
* Fraction of charged residues
* Fraction of positively minus negatively
charged residues.
* Score on how hydrophobic the surface is (-1 for every hydrophobic amino acid, normalised by protein size).
* Isoelectric point. 
* Charge of protein at pH 7. 
* Count of aromatic residues.
* The aromaticity of the protein
* Molecular weights in kDa. 
