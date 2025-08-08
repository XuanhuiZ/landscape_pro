from Bio.PDB import PDBParser
import numpy as np
import pandas as pd


def get_residue_coords(structure, chain_id, res_seq, icode=' '):
    """
    Extracts a (N_atoms × 3) array of atomic coords for one residue.
    
    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        Parsed PDB structure.
    chain_id : str
        Chain identifier (e.g. 'A').
    res_seq : int
        Residue sequence number.
    icode : str, optional
        Insertion code, if any (default=' ').
    
    Returns
    -------
    numpy.ndarray of shape (N_atoms, 3)
    """
    chain = structure[0][chain_id]
    # the key for a residue is a tuple: (hetfield, resseq, icode)
    res_id = (' ', res_seq, icode)
    residue = chain[res_id]
    coords = [atom.get_coord() for atom in residue.get_atoms()]
    return np.array(coords)

def min_distance_between(res1_coords, res2_coords):
    """
    Compute the minimum pairwise distance between two sets of coords.
    
    Parameters
    ----------
    res1_coords, res2_coords : numpy.ndarray
        Shape (N1,3) and (N2,3).
    
    Returns
    -------
    float
        Minimum Euclidean distance.
    """
    # broadcasting trick: (N1,1,3) - (1,N2,3) → (N1,N2,3)
    diffs = res1_coords[:, None, :] - res2_coords[None, :, :]
    d2 = np.sum(diffs*diffs, axis=2)
    return np.sqrt(d2.min())

def calculate_all_pairwise_min_dists(pdb_file, chain_id, res_list):
    """
    For a list of residue numbers, compute the min-distance for each pair.
    
    Parameters
    ----------
    pdb_file : str
        Path to .pdb file.
    chain_id : str
        Chain identifier.
    res_list : list of int
        Residue sequence numbers to compare.
    
    Returns
    -------
    dict of ((i,j) -> float)
        Keys are tuples of residue numbers, values are min distances.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('X', pdb_file)
    
    # preload all coords
    coords = {
        res: get_residue_coords(struct, chain_id, res)
        for res in res_list
    }
    
    dists = {}
    for i, r1 in enumerate(res_list):
        for r2 in res_list[i+1:]:
            d = min_distance_between(coords[r1], coords[r2])
            dists[(r1, r2)] = d
    return dists

def calculate_distance_matrix(pdb_file, chain_id, res_list):
    """
    Build an N×N matrix of min inter‐atomic distances between each pair of residues.

    Parameters
    ----------
    pdb_file : str
        Path to your .pdb file.
    chain_id : str
        Chain identifier (e.g. 'A').
    res_list : list of int
        Residue numbers (resseq) to include, in the order you want.

    Returns
    -------
    pandas.DataFrame
        Square DataFrame, indexed & columned by the residue numbers,
        whose (i,j) entry is the minimum atom–atom distance between
        res_list[i] and res_list[j] (in Å).
    """
    # parse once
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('X', pdb_file)

    # grab coords for each residue
    coords = {}
    for res in res_list:
        coords[res] = get_residue_coords(struct, chain_id, res)

    n = len(res_list)
    mat = np.zeros((n, n), dtype=float)

    # fill upper triangle (including diag)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                d = 0.0
            else:
                d = min_distance_between(coords[res_list[i]], coords[res_list[j]])
            mat[i, j] = d
            mat[j, i] = d  # mirror

    # wrap in DataFrame for nicer labeling
    return pd.DataFrame(mat,
                        index=[f"{chain_id}{r}" for r in res_list],
                        columns=[f"{chain_id}{r}" for r in res_list])


pdb_path = "/Users/zhouxuanhui/Desktop/landscape_pro/data/PTE.pdb"
chain = "A"
key_residues = [233, 254, 271, 272, 313, 306]

dist_df = calculate_distance_matrix(pdb_path, chain, key_residues)
dist_df.to_csv("/Users/zhouxuanhui/Desktop/landscape_pro/data/PTE_dis.csv")