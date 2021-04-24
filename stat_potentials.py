import numpy as np
import glob
import Bio.PDB as PDB

def matrix_retriever(structure): # , model, chain
    """Function that takes a structure as an input and outputs an matrix of the alfa carbon corrdinates"""
    coordinates_x = np.array(0)  # initialise empty arrays
    coordinates_y = np.array(0)
    coordinates_z = np.array(0)
    names = np.array(0)
    # parse the pdb file

    for model in structure:
        for chain in model:
            for residue in chain: #['A']: #[model][chain]:  # to take the residues from your structure, model 0, chain A which is the protein chain
                #if PDB.is_aa(residue.get_resname()):  # take out all the non amino acid residue
                #print(residue["CA"])
                try:
                    coordinates_x = np.append(coordinates_x, residue["CA"].get_coord()[0])  # get the coordinates
                    coordinates_y = np.append(coordinates_y, residue["CA"].get_coord()[1])  #
                    coordinates_z = np.append(coordinates_z, residue["CA"].get_coord()[2])
                    names = np.append(residue.get_resname())

                except:
                    pass

    coordinates_x = np.delete(coordinates_x, 0)  # delete first 0 in the matrix
    coordinates_y = np.delete(coordinates_y, 0)
    coordinates_z = np.delete(coordinates_z, 0)

    coordinates = np.column_stack([coordinates_x, coordinates_y, coordinates_z])  # merge the three matrices'''
    return coordinates, names

# may round then
# carbon

#coordinates_model0 = matrix_retriever("data/training/crystal_structs/A0A140NA.pdb", 0, "A")


def distance(a, b):
    """Function that calculates the distance between two coordinates"""
    dist = np.linalg.norm(a-b)


def statistical_potentials():
    """Function that calculates the statistical potential of a amino acid under sa specif structure feature"""
    # - boltz * temp * ln ()
    # boltzmans constant
    k = 5.67 * 10 ** (-8)
    # absolute temperature at room temperature
    T = 298.15  # K
    w = - k * T * np.log(prob_c_s / (prob_c * prob_s))
    return w


def call_statistical_potentials(filenames, aas):
    """Function that calculates probabilities before
        calculating the staistical potentials"""
    # 1 amino acid, 1 distance
    # find amino acids - calculate distance
    for file in filenames:
        # parse the pdb file
        p = PDB.PDBParser(QUIET=True)
        s = p.get_structure("", file)

        coordinates, names = matrix_retriever(s)
        #print(coordinates)

    return coordinates

amino_acids = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K',
                   'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

files = glob.glob("data/training/crystal_structs/*.pdb")
call_statistical_potentials(files, amino_acids)