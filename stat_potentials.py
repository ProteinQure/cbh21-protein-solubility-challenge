def statistical_potentials():
    """Function that calculates the statistical potential of a amino acid under sa specif structure feature"""
    # - boltz * temp * ln ()
    # boltzmans constant
    k = 5.67 * 10 ** (-8)
    # absolute temperature at room temperature
    T = 298.15  # K

    w = - k * T * np.log(prob_c_s / (prob_c * prob_s)

def call_statistical_potentials():
    """Function that calculates probabilities before
        calculating the staistical potentials"""
    # 1 amino acid, 1 distance
    # find amino acids - calculate distance


def matrix_retriever(structure, model, chain):
    """Function that takes a structure as an input and outputs an matrix of the alfa carbon corrdinates"""
    coordinates_x = np.array(0)  # initialise empty arrays
    coordinates_y = np.array(0)
    coordinates_z = np.array(0)
    names = np.array(0)

    for residue in structure[model][
        chain]:  # to take the residues from your structure, model 0, chain A which is the protein chain
        if PDB.is_aa(residue.get_resname()):  # take out all the non amino acid residue
            coordinates_x = np.append(coordinates_x, residue["CA"].get_coord()[0])  # get the coordinates
            coordinates_y = np.append(coordinates_y, residue["CA"].get_coord()[1])  #
            coordinates_z = np.append(coordinates_z, residue["CA"].get_coord()[2])
            names = residue.get_resname()

    coordinates_x = np.delete(coordinates_x, 0)  # delete first 0 in the matrix
    coordinates_y = np.delete(coordinates_y, 0)
    coordinates_z = np.delete(coordinates_z, 0)

    coordinates = matrix([coordinates_x, coordinates_y, coordinates_z, names])  # merge the three matrices
    return coordinates

# may round then
# carbon

#coordinates_model0 = matrix_retriever("data/training/crystal_structs/A0A140NA.pdb", 0, "A")


def distance(a, b):
    """Function that calculates the distance between two coordinates"""
    dist = np.linalg.norm(a-b)