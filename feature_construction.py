import freesasa
import Bio.PDB as PDB 
import Bio.PDB.DSSP as DSSP 
import glob
import numpy as np
from calc_fractions import *

def calc_length(filenames):

    prot_lengths = []

    for filename in filenames:
        # parse the pdb file
        p = PDB.PDBParser(QUIET=True)
        s = p.get_structure("", filename)

        # get the sequence length
        seq = 0
        for chain in s.get_chains():
            seq += len([_ for _ in chain.get_residues() if PDB.is_aa(_)])

        # save into numpy sheet rsa
        prot_lengths.append(seq)

    return prot_lengths
        

def calc_surface(filenames):

    surfaces = []
    for filename in filenames:

        # get the surface area
        structure = freesasa.Structure(filename)
        result = freesasa.calc(structure)
        area_classes = freesasa.classifyResults(result, structure)
        
        # save this into numpy sheet result.totalArea() 
        surface = result.totalArea()
        surfaces.append(surface)

    return surfaces


def ss_depth(filenames):

    frac_mod_beta_list = []
    frac_mod_alfa_list = []
    frac_exp_alfa_list = []

    for filename in filenames:
        p = PDB.PDBParser()
        structure = p.get_structure("", filename)
            
        model = structure[0]
        dssp = DSSP(model, filename)

        # DSSP data is accessed by a tuple (chain_id, res_id)
        all_residues = list(dssp.keys())

        # get a list with ASA values of all residues
        asa = [dssp[i][3] for i in all_residues]

        # 3 categories, like in paper
        burried = [0 if i <= 0.2 else 2 if i >= 0.5 else 1 for i in asa]

        # get a list with secondary structure of all residues (Q8)
        secondary_q8 = [dssp[i][2] for i in all_residues]

        # map the Q8 to Q3
        # helix = H, G, I
        # beta = B, E 
        # loop = rest 
        # 0 is alpha, 1 is beta, 2 is coil
        secondary_q3 = [0 if i in ['H', 'G', 'I'] else 1 if i in ['B', 'E'] else 2 for i in secondary_q8]

        # get the total number of sheets and helices
        count_helices = secondary_q3.count(0)
        count_sheets = secondary_q3.count(1)

        if count_sheets == 0:
            frac_mod_beta_list.append(0)
        else:
            # calculate fraction of moderately buried beta residues
            mod_beta = 0

            for i in range(len(burried)): 
                if burried[i] == 1 and secondary_q3[i] == 1:
                    mod_beta += 1
                
            frac_mod_beta = mod_beta / count_sheets
            frac_mod_beta_list.append(frac_mod_beta)

        if count_helices == 0:
            frac_mod_alfa_list.append(0)
            frac_exp_alfa_list.append(0)
        else:
            # calc fraction of moderately buried alfa residues
            mod_alfa = 0

            for i in range(len(burried)): 
                if burried[i] == 1 and secondary_q3[i] == 0:
                    mod_alfa += 1
                
            frac_mod_alfa = mod_alfa / count_helices
            frac_mod_alfa_list.append(frac_mod_alfa)

            # calc fraction of exposed a residues
            exp_alfa = 0

            for i in range(len(burried)): 
                if burried[i] == 2 and secondary_q3[i] == 0:
                    exp_alfa += 1
                
            frac_exp_alfa = exp_alfa / count_helices
            frac_exp_alfa_list.append(frac_exp_alfa)

    return frac_mod_beta_list, frac_mod_alfa_list, frac_exp_alfa_list





def compute_features(filenames):
    """"Takes list of pdb filenames as input and returns list of features"""

    ###### Protein IDs
    protIDs = [filename.split('/')[-1].split('.pdb')[0] for filename in filenames]
    
    ###### Protein Sequence Length
    print('Calculating Protein Sequence Length')
    prot_lengths = calc_length(filenames)

    ###### Surface Area
    print("Calculating Surface Area")
    surfaces = calc_surface(filenames)

    ###### Relative Surface Area
    print("Calculating Relative Surface Area")
    surface_seq = [surfaces[i] / prot_lengths[i] for i in range(len(surfaces))]

    ###### Secondary Structure
    print("Calculating Secondary Structure")
    frac_mod_beta_list, frac_mod_alfa_list, frac_exp_alfa_list = ss_depth(filenames)

    ###### Fractions of Negative and Positive
    print("Calculating Fractions of Negative and Positive")
    feats=feat_list(filenames)
    # turn list of tuples into tuple of lists
    frac_k_minus_r, frac_neg, frac_pos, frac_charged, pos_minus_neg= ([a for a,b,c,d,e in feats],
                                                                      [b for a,b,c,d,e in feats],
                                                                      [c for a,b,c,d,e in feats],
                                                                      [d for a,b,c,d,e in feats],
                                                                      [e for a,b,c,d,e in feats])

    # Save features in file to pass to R script
    print("Saving features")
    arr = np.column_stack((protIDs, surfaces, prot_lengths, surface_seq, frac_mod_beta_list, frac_mod_alfa_list,
                           frac_exp_alfa_list,frac_k_minus_r,frac_neg,frac_pos,frac_charged,pos_minus_neg))

    print(arr)
    np.savetxt("features.csv", arr, delimiter=",")


if __name__ == "__main__":
    filenames = glob.glob("data/training/crystal_structs/*.pdb")

    compute_features(filenames)