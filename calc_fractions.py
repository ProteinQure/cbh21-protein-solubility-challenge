import Bio.PDB.DSSP as DSSP
import Bio.PDB as PDB
import numpy as np
import glob
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

def get_feats(file):

    alphas=['H','I','G']
    betas=['B','E']

    p = PDB.PDBParser(QUIET=True)

    dssp_tuple = dssp_dict_from_pdb_file(file)
    dssp_dict = dssp_tuple[0]

    all_residues = list(dssp_dict.keys())

    # get total fraction of alpha and beta
    alph_frac=0
    bet_frac=0
    for amino_acid in all_residues:
        ss = dssp_dict[amino_acid][1]
        if ss in alphas:
            alph_frac+=1
        if ss in betas:
            bet_frac+=1


    len_prot=len(dssp_dict)

    # looking at beta residues:
    tot_bet_bur=0
    tot_bet_mod = 0
    for amino_acid in all_residues:
        ss=dssp_dict[amino_acid][1]
        RASA=dssp_dict[amino_acid][2]
        if ss in betas:
            if RASA<100:
                tot_bet_bur+=1
            elif RASA<150:
                tot_bet_mod += 1

    # calculate fraction of buried beta residues, append to list
    try:
        bet_bur=tot_bet_bur/bet_frac
    except:
        bet_bur=0
    # calculate fraction of moderately buried beta residues, append to list
    try:
        bet_mod=(tot_bet_mod/bet_frac)
    except:
        bet_mod=0

    # looking at alpha residues:
    tot_al_mod = 0
    tot_al_exp = 0
    for amino_acid in all_residues:
        ss = dssp_dict[amino_acid][1]
        RASA = dssp_dict[amino_acid][2]
        if ss in alphas:
            if RASA > 150:
                tot_al_exp += 1
            elif 100 < RASA < 150:
                tot_al_mod += 1
        # calculate fraction of moderately buried alpha residues, append to list
    try:
        al_mod=(tot_al_mod / alph_frac)
    except:
        al_mod=0
    # calculate fraction of moderately buried beta residues, append to list
    try:
        al_exp=(tot_al_exp / alph_frac)
    except:
        al_exp='NA'


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
    try:
        frac_k_minus_r=aas['K']-aas['R']
    except:
        try:
            frac_k_minus_r = -aas['R']
        except:
            try:
                frac_k_minus_r = aas['K']
            except:
                frac_k_minus_r = 0

    # fraction of negatively charged residues
    negs=['D','E']
    frac_neg=0
    for neg in negs:
        try:
            frac_neg+=aas[neg]/len_prot
        except:
            pass

    # fraction of positively charged residues
    poss=['K','H','R']
    frac_pos=0
    for pos in poss:
        try:
            frac_pos+=aas[pos]/len_prot
        except:
            pass

    # fraction of charged residues
    charged=['D','E','K','H','R']
    frac_charged = 0
    for ch in charged:
        try:
            frac_charged += aas[ch] / len_prot
        except:
            pass

    # fraction of positively minus negatively charged residues
    pos_minus_neg=frac_pos-frac_neg

    fracs=(#bet_bur,
                #bet_mod,
                #al_mod,
                #al_exp,
                #aas,
                frac_k_minus_r,
                frac_neg,
                frac_pos,
                frac_charged,
                pos_minus_neg)

    return(fracs)

def feat_list(list_of_proteins):
    list_feats=[]
    for file in list_of_proteins:
        feats=get_feats(file)
        list_feats.append(feats)
    return(list_feats)

# list_prots=glob.glob("data/training/crystal_structs/*.pdb")
