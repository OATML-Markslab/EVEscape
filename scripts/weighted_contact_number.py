import pandas as pd
from Bio.Data import SCOPData
from Bio.PDB import PDBParser
from Bio.PDB import is_aa
import warnings

# Adapted from ProteinER (See acknowledgements for details/license info)


def inv_sq_distance(coord1, coord2):
    '''
    Returns the inverse squared distance between any two coordinates.
    '''
    distance = 0.0
    for i, j in zip(coord1, coord2):
        distance += (i - j)**2
    return 1 / distance


def calculate_wcn(residues):
    '''
    Calculates weighted contact number (WCN).
    '''
    for residue in residues:
        wcn_ca = 0
        wcn_sc = 0
        for other_residue in residues:
            if residue != other_residue:
                wcn_ca += inv_sq_distance(residue['coord_ca'],
                                          other_residue['coord_ca'])
                wcn_sc += inv_sq_distance(residue['sidechain_center'],
                                          other_residue['sidechain_center'])
        residue['wcn_ca'] = wcn_ca
        residue['wcn_sc'] = wcn_sc

    return residues


def process_residue(residue):
    '''
    Processes a single residue to determine the coordinates of the alpha-carbon
    and the sidechain center-of-mass. Also checks for missing atoms in a
    residue.
    '''
    output_dict = {}
    atoms_seen = []
    # Convert three letter amino acid to one letter
    output_dict['pdb_aa'] = SCOPData.protein_letters_3to1[residue.resname]
    # Grab residue number AND any insertion site labeling (11A, 11B, etc.)
    output_dict['pdb_position'] = str(residue.get_id()[1]) + \
                                  residue.get_id()[2].strip()
    output_dict['chain'] = residue.get_full_id()[2]
    # Coordinates of all sidechain atoms in this residue
    sidechain_coords = []
    for atom in residue:
        atoms_seen.append(atom.name)
        if atom.name == 'CA':
            # Save alpha-carbon coordinates
            output_dict['coord_ca'] = atom.get_coord()
        if atom.name not in ['C', 'CA', 'O', 'N']:
            # Must be a sidechain atom...
            sidechain_coords.append(atom.get_coord())

    warning_message = "Missing {} in residue (" + \
                        str(output_dict['pdb_position']) + ", " + \
                        str(output_dict['pdb_aa']) + ")"

    for mainchain_atom in ['N', 'C', 'O']:
        # Warn about any missing mainchain atoms
        if mainchain_atom not in atoms_seen:
            warnings.warn(warning_message.format(mainchain_atom),
                          RuntimeWarning)
    if 'coord_ca' not in output_dict:
        # Cannot calculate WCN without at least alpha-carbon
        raise RuntimeError(
            warning_message.format('CA') + '. Cannot calculate C-alpha WCN.')

    if len(sidechain_coords) == 0:
        # Warn about missing sidechain for amino acids other than glycine
        if output_dict['pdb_aa'] != 'G':
            warnings.warn(
                warning_message.format('sidechain') + '. Using CA instead.',
                RuntimeWarning)
        sidechain_coords.append(output_dict['coord_ca'])

    # Calculate side chain center of mass
    output_dict['sidechain_center'] = sum(sidechain_coords)/\
                                      len(sidechain_coords)

    return output_dict


def collect_coordinates(structure):
    '''
    Loops over all residues in a structure and collects coordinates for alpha-
    carbons and sidechain center-of-mass. Returns a list of dictionaries, where
    each dictionary corresponds to residue in the structure.
    '''
    output_list = []
    for residue in structure.get_residues():
        if is_aa(residue):
            output_list.append(process_residue(residue))
    return output_list


def create_wcn_file(pdb_path, chain_id):
    # Load in PDB with BioPython
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure("pdb", pdb_path)

    # Collect coordinate information
    output_list = collect_coordinates(structure)

    output_list = [o for o in output_list if o["chain"] in chain_id]

    # Calculate WCN from coordinates
    output_list = calculate_wcn(output_list)
    return (pd.DataFrame(output_list))


def add_wcn_to_site_annotations(pdb_path, chain_id):
    wcns = create_wcn_file(pdb_path, chain_id)

    wcns = wcns[["chain", "pdb_aa", "pdb_position", "wcn_ca", "wcn_sc"]]

    return (wcns)
