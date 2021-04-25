import Bio
import os
import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import freesasa


def get_area_classes(file):
    struct = freesasa.Structure(file)
    result = freesasa.calc(struct)
    area_classes = freesasa.classifyResults(result,struct)
    list_areas = [(list(area_classes.values())[0]),
                  (list(area_classes.values())[1]),
                  result.totalArea()]
    return list_areas
    #dict = {"polar":"value", "Apolar":"Value"}

path  = r'C:\Users\Kai Armstrong\PycharmProjects\cbh21-protein-solubility-challenge-main\cbh21-protein-solubility-challenge-main\data\test\yeast_crystal_structs'


#solubility_values_df = pd.read_csv('solubility_values.csv')

parser=PDBParser()

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

#struct = freesasa.Structure('P0ADV9.pdb')
#result = freesasa.calc(struct)
#area_classes = freesasa.classifyResults(result, struct)
#print(result.totalArea())
#print(area_classes)



for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            areas = get_area_classes(entry)
            polar_area.append(areas[0])
            apolar_area.append(areas[1])
            total_area.append(areas[2])

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            for record in SeqIO.parse(entry, "pdb-atom"):
                sequence = str(record.seq).replace('X', 'G')
                protein = ProteinAnalysis(sequence)
                p_len.append(len(sequence))
                mol_w.append(protein.molecular_weight())
                iso_p.append(protein.isoelectric_point())
                smell.append(protein.aromaticity())
                taste_factor.append(protein.gravy())
                insta_ind.append(protein.instability_index())
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

#name_col = solubility_values_df['protein'].tolist()
#sol_col = solubility_values_df['solubility'].tolist()

#list_lists.append(name_col)
for values_count in d_count.values():
    list_lists.append(list(values_count))
for values_count in d_perc.values():
    list_lists.append(list(values_count))
list_lists.append(p_len)
list_lists.append(smell)
list_lists.append(taste_factor)
list_lists.append(mol_w)
list_lists.append(iso_p)
list_lists.append(insta_ind)
list_lists.append(helter_skeler)
list_lists.append(turnip)
list_lists.append(garfield)
list_lists.append(polar_area)
list_lists.append(apolar_area)
list_lists.append(total_area)
#list_lists.append(sol_col)

df = pd.DataFrame(list_lists)
df = df.transpose()
df.to_csv('yeast_crystal.csv', index=False)


"""
def get_sequence(record):
    sequence = str(record.seq)
    return sequence

def get_prot(sequence):
    protein = ProteinAnalysis(sequence)
    return protein

def get_molw(protein):
     mol_w.append(protein.molecular_weight())
filenames = []
ids = []
counter = 0
for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            for record in SeqIO.parse(entry, "pdb-atom"):
                try:
                    get_molw(get_prot(get_sequence(record)))
                except:
                    break
                counter += 1
                ids.append(record.id)
                print(mol_w)
print(counter)
for entry1 in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry1)):
        if entry1.endswith('.pdb'):
            filenames.append(entry1)
print(set(filenames)-set(ids))
#pdbfile=open('A0A140NA.pdb')
#p = PDBParser()
#structure = p.get_structure("A0A140NA", "cbh21-protein-solubility-challenge-main/data/training/crystal_structs/A0A140NA.pdb")
#model = structure[0]
#dssp = DSSP(model, "/cbh21-protein-solubility-challenge-main/data/training/crystal_structs/A0A140NA.pdb")


pdbfile=open('A0A140NA.pdb')
p = PDBParser()
s=p.get_structure('X',pdbfile)
m=s[0]
RADIUS=12.0
hse=HSExposureCB(m, RADIUS)
residue_list=Selection.unfold_entities(m,'R')
hse=HSExposureCA(m, RADIUS)
residue_list = Selection.unfold_entities(m,'R')
for r in residue_list[:10]:
    print(r.get_resname(), r.xtra)
for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            areas = get_area_classes(entry)
            polar_area.append(areas[0])
            apolar_area.append(areas[1])
            total_area.append(areas[2])

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            for record in SeqIO.parse(entry, "pdb-atom"):
                sequence = str(record.seq).replace('X', 'G')
                protein = ProteinAnalysis(sequence)
                p_len.append(len(sequence))
                mol_w.append(protein.molecular_weight())
                iso_p.append(protein.isoelectric_point())
                smell.append(protein.aromaticity())
                taste_factor.append(protein.gravy())
                insta_ind.append(protein.instability_index())
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

name_col = solubility_values_df['protein'].tolist()
sol_col = solubility_values_df['solubility'].tolist()

list_lists.append(name_col)
for values_count in d_count.values():
    list_lists.append(list(values_count))
for values_count in d_perc.values():
    list_lists.append(list(values_count))
list_lists.append(p_len)
list_lists.append(smell)
list_lists.append(taste_factor)
list_lists.append(mol_w)
list_lists.append(iso_p)
list_lists.append(insta_ind)
list_lists.append(helter_skeler)
list_lists.append(turnip)
list_lists.append(garfield)
list_lists.append(polar_area)
list_lists.append(apolar_area)
list_lists.append(total_area)
list_lists.append(sol_col)

df = pd.DataFrame(list_lists)
df = df.transpose()
df.to_csv('data_features_train.csv', index=False)

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            areas = get_area_classes(entry)
            polar_area.append(areas[0])
            apolar_area.append(areas[1])
            total_area.append(areas[2])

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            for record in SeqIO.parse(entry, "pdb-atom"):
                sequence = str(record.seq).replace('X', 'G')
                protein = ProteinAnalysis(sequence)
                p_len.append(len(sequence))
                mol_w.append(protein.molecular_weight())
                iso_p.append(protein.isoelectric_point())
                smell.append(protein.aromaticity())
                taste_factor.append(protein.gravy())
                insta_ind.append(protein.instability_index())
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

name_col = solubility_values_df['protein'].tolist()
sol_col = solubility_values_df['solubility'].tolist()

list_lists.append(name_col)
for values_count in d_count.values():
    list_lists.append(list(values_count))
for values_count in d_perc.values():
    list_lists.append(list(values_count))
list_lists.append(p_len)
list_lists.append(smell)
list_lists.append(taste_factor)
list_lists.append(mol_w)
list_lists.append(iso_p)
list_lists.append(insta_ind)
list_lists.append(helter_skeler)
list_lists.append(turnip)
list_lists.append(garfield)
list_lists.append(polar_area)
list_lists.append(apolar_area)
list_lists.append(total_area)
list_lists.append(sol_col)

df = pd.DataFrame(list_lists)
df = df.transpose()
df.to_csv('data_features_train.csv', index=False)

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            areas = get_area_classes(entry)
            polar_area.append(areas[0])
            apolar_area.append(areas[1])
            total_area.append(areas[2])

for entry in os.listdir(path):
    if os.path.isfile(os.path.join(path, entry)):
        if entry.endswith('.pdb'):
            for record in SeqIO.parse(entry, "pdb-atom"):
                sequence = str(record.seq).replace('X', 'G')
                protein = ProteinAnalysis(sequence)
                p_len.append(len(sequence))
                mol_w.append(protein.molecular_weight())
                iso_p.append(protein.isoelectric_point())
                smell.append(protein.aromaticity())
                taste_factor.append(protein.gravy())
                insta_ind.append(protein.instability_index())
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

name_col = solubility_values_df['protein'].tolist()
sol_col = solubility_values_df['solubility'].tolist()

list_lists.append(name_col)
for values_count in d_count.values():
    list_lists.append(list(values_count))
for values_count in d_perc.values():
    list_lists.append(list(values_count))
list_lists.append(p_len)
list_lists.append(smell)
list_lists.append(taste_factor)
list_lists.append(mol_w)
list_lists.append(iso_p)
list_lists.append(insta_ind)
list_lists.append(helter_skeler)
list_lists.append(turnip)
list_lists.append(garfield)
list_lists.append(polar_area)
list_lists.append(apolar_area)
list_lists.append(total_area)
list_lists.append(sol_col)

df = pd.DataFrame(list_lists)
df = df.transpose()
df.to_csv('data_features_train.csv', index=False)
"""
