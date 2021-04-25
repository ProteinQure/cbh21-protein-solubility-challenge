import Bio
import os
import numpy as np
import pandas as pd
import temppathlib
import zipfile
import argparse
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import freesasa
from Bio import PDB
import pickle


def get_area_classes(file):
    struct = freesasa.Structure(file)
    result = freesasa.calc(struct)
    area_classes = freesasa.classifyResults(result,struct)
    list_areas = [(list(area_classes.values())[0]),
                  (list(area_classes.values())[1]),
                  result.totalArea()]
    return list_areas
    #dict = {"polar":"value", "Apolar":"Value"}
"""
def get_seq(pdbfile):
    file = str(pdbfile)
    p = PDBParser(PERMISSIVE=0)
    structure = p.get_structure(file[-3], pdbfile)
    ppb = PPBuilder()
    seq = ''
    for pp in ppb.build_peptides(structure):
        seq += pp.get_sequence()
    return seq
"""
#parser=PDBParser()

amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

list_lists = []

d_count = {}
for i in amino_acids:
    d_count["count_{0}".format(i)] = []
d_perc = {}
for j in amino_acids:
    d_perc["perc_{0}".format(j)] = []

p_len = [] #length of protein
smell = [] #aromaticity
taste_factor = [] # < 0 = hydrophilic, > 0 = hydrophobic -> gravy
mol_w = [] #molecular weight
iso_p = [] #isoelectric point
insta_ind = [] #under 40 = stable (ish)
helter_skeler = [] #helices
turnip = [] #turns
garfield = [] #sheets
polar_area = [] #polar surface area
apolar_area = [] #apolar surface area
total_area = [] #total surface area
char_at_acid = []
char_at_neutral = []
char_at_base = []

parser = argparse.ArgumentParser()
parser.add_argument("--infile", type=str, default="data/test.zip")
parser.add_argument("--model", type=str, default="model.pkl")
args = parser.parse_args()

#protein_parser = PDBParser()

with temppathlib.TemporaryDirectory() as tmpdir:
    # unzip the file with all the test PDBs
    with zipfile.ZipFile(args.infile, "r") as zip_:
        zip_.extractall(tmpdir.path)
        for test_pdb in tmpdir.path.glob("*.pdb"):
            struct = freesasa.Structure(str(test_pdb))
            result = freesasa.calc(struct)
            areas_classes = freesasa.classifyResults(result,struct)
            list_areas = [(list(areas_classes.values())[0]),
                          (list(areas_classes.values())[1]),
                          result.totalArea()]

            polar_area.append(list_areas[0])
            apolar_area.append(list_areas[1])
            total_area.append(list_areas[2])

print('done')
with temppathlib.TemporaryDirectory() as tmpdir:
        # unzip the file with all the test PDBs
        with zipfile.ZipFile(args.infile, "r") as zip_:
            zip_.extractall(tmpdir.path)

            for test_pdb in tmpdir.path.glob("*.pdb"):
                for record in SeqIO.parse(test_pdb, "pdb-atom"):
                    sequence = str(record.seq).replace('X', 'G')
                    protein = ProteinAnalysis(str(sequence))
                    p_len.append(len(sequence))
                    mol_w.append(protein.molecular_weight())
                    iso_p.append(protein.isoelectric_point())
                    smell.append(protein.aromaticity())
                    taste_factor.append(protein.gravy())
                    insta_ind.append(protein.instability_index())
                    char_at_acid.append(protein.charge_at_pH(1))
                    char_at_neutral.append(protein.charge_at_pH(7))
                    char_at_base.append(protein.charge_at_pH(14))
                    helter_skeler.append(protein.secondary_structure_fraction()[0])
                    turnip.append(protein.secondary_structure_fraction()[1])
                    garfield.append(protein.secondary_structure_fraction()[2])
                    for x in amino_acids:
                        n = protein.count_amino_acids()[x]
                        for y in d_count.keys():
                            if y[-1] == x:
                                d_count[y].append(n)
                    for a in amino_acids:
                        m = protein.get_amino_acids_percent()[a]
                        for b in d_perc.keys():
                            if b[-1] == a:
                                d_perc[b].append(m)
                #areas = get_area_classes(test_pdb)
                #polar_area.append(areas[0])
                #apolar_area.append(areas[1])
                #total_area.append(areas[2])
print('done')

for values_count in d_count.values():
    list_lists.append(list(values_count))
for values_count in d_perc.values():
    list_lists.append(list(values_count))
list_lists.append(p_len)
list_lists.append(smell)
list_lists.append(taste_factor)
list_lists.append(mol_w)
list_lists.append(iso_p)
list_lists.append(char_at_acid)
list_lists.append(char_at_neutral)
list_lists.append(char_at_base)
list_lists.append(insta_ind)
list_lists.append(helter_skeler)
list_lists.append(turnip)
list_lists.append(garfield)
list_lists.append(polar_area)
list_lists.append(apolar_area)
list_lists.append(total_area)

df = pd.DataFrame(list_lists)
df = df.transpose()
df.reset_index(drop=True, inplace=True)
print(df)

df.to_csv('features_from_zip.csv', index=False)

continuous_data = df.columns.tolist()

print(continuous_data)

x_test = pd.read_csv('features_from_zip.csv')
print(x_test)

protein_names = []

with temppathlib.TemporaryDirectory() as tmpdir:
        # unzip the file with all the test PDBs
        with zipfile.ZipFile(args.infile, "r") as zip_:
            zip_.extractall(tmpdir.path)
            for test_pdb in tmpdir.path.glob("*.pdb"):
              filename = os.path.basename(test_pdb)
              protein_names.append(str(filename)[:-4])

print(protein_names)

with open(args.model, 'rb') as f:
    estimator = pickle.load(f)
    print(estimator)
    predictions = estimator.predict(x_test)
    print(predictions)
    y_pred = pd.DataFrame(list(zip(protein_names,predictions)), columns=['protein','solubility']).to_csv('predictions.csv', index=False, header=True)