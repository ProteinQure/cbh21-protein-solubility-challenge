import Bio.PDB.DSSP as DSSP
import Bio.PDB as PDB
import numpy as np
import glob
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

def get_feats(file):

    p = PDB.PDBParser(QUIET=True)
    structure = p.get_structure(file, file)

    dssp_tuple = dssp_dict_from_pdb_file(file)
    dssp_dict = dssp_tuple[0]


    #dssp = DSSP(model, "data/training/crystal_structs/A0A140NA.pdb")
    # DSSP data is accessed by a tuple (chain_id, res_id)
    #a_key = list(dssp_dict.keys())[2]

    all_residues = list(dssp_dict.keys())

    dssp_info = [dssp_dict[i] for i in all_residues]

    len_prot=len(dssp_info)

    # looking at beta residues:
    tot_bet_bur=0
    tot_bet_mod = 0
    for amino_acid in all_residues:
        ss=dssp_dict[amino_acid][1]
        RASA=dssp_dict[amino_acid][2]
        if ss=='B' or ss=='E':
            if RASA<100:
                tot_bet_bur+=1
            elif RASA<150:
                tot_bet_mod += 1

    # calculate fraction of buried beta residues, append to list
    bet_bur=tot_bet_bur/len_prot
    # calculate fraction of moderately buried beta residues, append to list
    bet_mod=(tot_bet_mod / len_prot)

    # looking at alpha residues:
    tot_al_mod = 0
    tot_al_exp = 0
    for amino_acid in all_residues:
        ss = dssp_dict[amino_acid][1]
        RASA = dssp_dict[amino_acid][2]
        if ss == 'H':
            if RASA > 150:
                tot_al_exp += 1
            elif 100 < RASA < 150:
                tot_al_mod += 1
        # calculate fraction of moderately buried alpha residues, append to list
        al_mod=(tot_al_mod / len_prot)
        # calculate fraction of moderately buried beta residues, append to list
        al_exp=(tot_al_exp / len_prot)


    # calc fraction of each of the 20 amino acid types
    aas={}
    for amino_acid in all_residues:
        aa = dssp_dict[amino_acid][0]
        if aa in aas.keys():
            aas[aa]+=1
        else:
            aas[aa]=0
    for aa in aas.keys():
        aas[aa]=aas[aa]/len_prot

    # calc fraction of K minus fraction of R
    frac_k_minus_r=aas['K']-aas['R']

    # fraction of negatively charged residues


    # fraction of charged residues

    # fraction of positively minus negatively charged residues

    list_fracs=[bet_bur,
                      bet_mod,
                      al_mod,
                      al_exp,
                aas,
                frac_k_minus_r]
    print(list_fracs)
    return(list_fracs)

for file in glob.glob("data/training/crystal_structs/*.pdb"):
    get_feats(file)