import Bio.PDB.DSSP as DSSP
import Bio.PDB as PDB
import numpy as np
import glob
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


p = PDB.PDBParser()

def get_feats(file):
    p = PDB.PDBParser(QUIET=True)
    structure = p.get_structure(file, file)

    model = structure[0]

    dssp_tuple = dssp_dict_from_pdb_file(file)
    dssp_dict = dssp_tuple[0]


    #dssp = DSSP(model, "data/training/crystal_structs/A0A140NA.pdb")
    # DSSP data is accessed by a tuple (chain_id, res_id)
    #a_key = list(dssp_dict.keys())[2]

    all_residues = list(dssp_dict.keys())

    dssp_info = [dssp_dict[i] for i in all_residues]

    len_prot=len(dssp_info)

    print(len_prot)
    tot_bet=0
    for amino_acid in all_residues:
        ss=dssp_dict[amino_acid][1]
        if ss=='B' or ss=='E':
            tot_bet+=1
    frac_bet=tot_bet/len_prot
    print(frac_bet)







    # calculate fraction of buried beta residues

    # calculate fraction of moderately buried beta residues

    # calc fraction of moderately buried alfa residues

    # calc fraction of exposed a residues

    # calc fraction of each of the 20 amino acid types

    # calc fraction of K minus fraction of R

    # fraction of negatively charged residues

    # fraction of charged residues

    # fraction of positively minus negatively charged residues

for file in glob.glob("data/training/crystal_structs/*.pdb"):
    get_feats(file)