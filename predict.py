"""
The entry point for your prediction algorithm.
"""

from __future__ import annotations
import argparse
import csv
import itertools
from pathlib import Path
import pprint
from typing import Any
import zipfile

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Structure import Structure
import temppathlib


def predict(pdb_file: Path) -> float:
    """
    The function that puts it all together: parsing the PDB file, generating
    features from it and performing inference with the ML model.
    """

    # parse PDB
    parser = PDBParser()
    structure = parser.get_structure(pdb_file.stem, pdb_file)

    # featurize + perform inference
    features = featurize(structure)
    predicted_solubility = ml_inference(features)

    return predicted_solubility


def featurize(structure: Structure) -> list[Any]:
    """
    Calculates 3D ML features from the `structure`.
    """

    # get all the residues
    residues = [res for res in structure.get_residues()]

    # calculate some random 3D features (you should be smarter here!)
    protein_length = residues[1]["CA"] - residues[-2]["CA"]
    angle = calc_dihedral(
        residues[1]["CA"].get_vector(),
        residues[2]["CA"].get_vector(),
        residues[-3]["CA"].get_vector(),
        residues[-2]["CA"].get_vector(),
    )
    # create the feature vector
    features = [protein_length, angle]

    return features


def ml_inference(features: list[Any]) -> float:
    """
    This would be a function where you normalize/standardize your features and
    then feed them to your trained ML model (which you would load from a file).
    """

    # this is my stupid manual ML model
    if features[0] > 15.0 and features[1] > 0.5:
        return 60
    elif features[0] > 30.0 and features[1] > 1.5:
        return 80

    return 20


if __name__ == "__main__":

    # set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, default="data/test.zip")
    args = parser.parse_args()

    predictions = []
    # use a temporary directory so we don't pollute our repo
    with temppathlib.TemporaryDirectory() as tmpdir:
        # unzip the file with all the test PDBs
        with zipfile.ZipFile(args.infile, "r") as zip_:
            zip_.extractall(tmpdir.path)

        # iterate over all test PDBs and generate predictions
        for test_pdb in tmpdir.path.glob("*.pdb"):
            predictions.append({"protein": test_pdb.stem, "solubility": predict(test_pdb)})

    # save to csv file, this will be used for benchmarking
    outpath = "predictions.csv"
    with open(outpath, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=["protein", "solubility"])
        writer.writeheader()
        writer.writerows(predictions)

    # print predictions to screen
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(predictions)
