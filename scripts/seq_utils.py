import numpy as np
from io import StringIO
import os, io, random, glob
import string
import pandas as pd
from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Data import SCOPData
from Bio.PDB import PDBParser
from Bio.PDB import is_aa

alphabet = "ACDEFGHIKLMNPQRSTVWY"

Codon_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "stop",
    "TAG": "stop",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "TGA": "stop",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

AminoAcid = list(set(Codon_AA.values()))


def fa_to_pd(filename):
    '''
    Read fasta file and return pandas dataframe with each row containing 
    a sequence header and sequence
    '''
    column_labels = ['header', 'Sequence']
    ali = []

    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    for fasta in fasta_sequences:
        seq_id, seq = fasta.id, str(fasta.seq)
        ali.append([seq_id, seq])

    ali = pd.DataFrame(ali, columns=column_labels)
    return ali


def write_fa(fa_table, outfile, headname, seqname):
    '''write seq table to fasta'''
    with open(outfile, 'w') as f:
        for i, row in fa_table.iterrows():
            f.write('>' + str(row[headname]) + '\n')
            f.write(row[seqname] + '\n')


def alphanumeric_index_to_numeric_index(
    n,
    alpha_to_numeric={
        a: (i + 1) / 100
        for i, a in enumerate(string.ascii_uppercase)
    }):
    '''convert pdb alphanumeric residue index to numeric index'''
    n = (''.join(c for c in n if c.isdigit())
         or None, ''.join(c for c in n if c.isalpha()) or None)

    if n[1]:
        return (int(n[0]) + alpha_to_numeric[n[1].upper()])
    else:
        return int(n[0])


def pairwise_align(seq1, seq2):
    '''Pairwise align protein sequences using biopython'''
    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.localds(seq1, seq2, blosum62, -10, -0.5)
    return (alignments[0][0], alignments[0][1])


def map_indices(seq_i, start_i, end_i, seq_j, start_j, end_j, gaps=("-", ".")):
    """
    Compute index mapping between positions in two
    aligned sequences
    Taken from evcouplings (https://github.com/debbiemarkslab/EVcouplings)
    Parameters
    ----------
    seq_i : str
        First aligned sequence
    start_i : int
        Index of first position in first sequence
    end_i : int
        Index of last position in first sequence
        (used for verification purposes only)
    seq_j : str
        Second aligned sequence
    start_j : int
        Index of first position in second sequence
    end_j : int
        Index of last position in second sequence
        (used for verification purposes only)
    Returns
    -------
    pandas.DataFrame
        Mapping table containing assignment of
        1. index in first sequence (i)
        2. symbol in first sequence (A_i)
        3. index in second sequence (j)
        4. symbol in second sequence (A_j)
    """
    NA = np.nan
    pos_i = start_i
    pos_j = start_j
    mapping = []

    for i, (res_i, res_j) in enumerate(zip(seq_i, seq_j)):
        # Do we match two residues, or residue and a gap?
        # if matching two residues, store 1 to 1 mapping.
        # Store positions as strings, since pandas cannot
        # handle nan values in integer columns
        if res_i not in gaps and res_j not in gaps:
            mapping.append([str(pos_i), res_i, str(pos_j), res_j])
        elif res_i not in gaps:
            mapping.append([str(pos_i), res_i, NA, NA])
        elif res_j not in gaps:
            mapping.append([NA, NA, str(pos_j), res_j])

        # adjust position in sequences if we saw a residue
        if res_i not in gaps:
            pos_i += 1

        if res_j not in gaps:
            pos_j += 1

    assert pos_i - 1 == end_i and pos_j - 1 == end_j

    return pd.DataFrame(mapping, columns=["i", "A_i", "j", "A_j"])


def remap_to_target_seq(new_seq, target_seq):
    '''
    Make mapping table with residue IDs and indices for a new sequence 
    relative to a target sequence
    '''
    aligned = pairwise_align(new_seq, target_seq)
    map_table = map_indices(aligned[0], 1, len(new_seq), aligned[1], 1,
                            len(target_seq))
    return (map_table)


def remap_pdb_seq_to_target_seq(pdb_file,
                                chain_list,
                                target_seq_file,
                                alphabet='ACDEFGHIKLMNPQRSTVWY'):
    '''Remap chains from a pdb file to a target sequence'''
    target_seq = fa_to_pd(target_seq_file).Sequence.values[0]
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure("pdb", pdb_file)
    output_list = []

    for i, residue in enumerate(structure.get_residues()):

        if is_aa(residue):
            output_dict = {}
            # Convert three letter amino acid to one letter
            output_dict['pdb_res'] = SCOPData.protein_letters_3to1[
                residue.resname]
            # Grab residue number AND any insertion site labeling (11A, 11B, etc.)
            output_dict['pdb_position'] = str(residue.get_id()[1]) + \
                                      residue.get_id()[2].strip()
            output_dict['chain'] = residue.get_full_id()[2]
            output_dict['ind'] = i
            output_list.append(output_dict)

    resis = pd.DataFrame(output_list)
    resis = resis[resis.chain.isin(chain_list)]
    resis['pdb_i'] = resis.sort_values(['chain', 'ind'
                                        ]).groupby('chain').cumcount() + 1
    map_tables = []

    for chain in chain_list:
        chain_seq = ''.join(resis[resis.chain == chain].sort_values(
            ['ind']).pdb_res.tolist())
        map_table = remap_to_target_seq(chain_seq, target_seq)
        map_table = map_table.rename(columns={
            'i': 'pdb_i',
            'A_i': 'pdb_res',
            'j': 'target_i',
            'A_j': 'target_res'
        })
        map_table['chain'] = chain
        map_table = map_table[~map_table.pdb_i.isna()]
        map_table = map_table[~map_table.target_i.isna()]
        map_table['pdb_i'] = map_table.pdb_i.astype(int)
        map_table['target_i'] = map_table.target_i.astype(int)
        map_tables.append(map_table)

    map_table = pd.concat(map_tables)
    resis = resis.drop(columns=['ind'])
    map_table = map_table.merge(resis,
                                how='left',
                                on=['chain', 'pdb_i', 'pdb_res'])
    return (map_table)


def make_mut_table(seqfile, alphabet=alphabet):
    '''Make DataFrame of all single mutations for a given protein sequence'''
    samples = []
    seq = fa_to_pd(seqfile).Sequence[0]

    for i, old in enumerate(seq):
        for new in alphabet:
            if new.upper() != old.upper():
                samples.append({
                    'i': (i + 1),
                    'wt': old.upper(),
                    'mut': new.upper()
                })

    return pd.DataFrame(samples)


def remap_struct_df_to_target_seq(struct_df, chainlist, map_table):
    '''
    remap dataframe derived from protein structure (like DSSP) to target 
    sequence, given mapping table
    '''
    struct_df = struct_df[struct_df['chain'].isin(chainlist)]
    struct_df['pdb_i'] = struct_df.sort_values(
        ['chain', 'i']).groupby('chain').cumcount() + 1
    struct_df = struct_df.rename(columns={
        'i': 'pdb_numbering',
        'wt': 'pdb_res'
    })
    struct_df = struct_df.merge(map_table,
                                how='left',
                                on=['chain', 'pdb_i', 'pdb_res'])
    struct_df = struct_df.drop(columns='pdb_i')
    struct_df = struct_df.rename(columns={'target_res': 'wt', 'target_i': 'i'})

    return (struct_df)


def find_syn_codons(a, Codon_AA):
    syn_cdns = []
    for c in Codon_AA:
        if Codon_AA[c] == a:
            syn_cdns.append(c)
    return (syn_cdns)


def syn_cdn_dict(AminoAcid, Codon_AA):
    AA_Codon = dict.fromkeys(AminoAcid)
    for a in AminoAcid:
        AA_Codon[a] = find_syn_codons(a, Codon_AA)
    return (AA_Codon)


def nuc_diff(source, target):
    '''
        Returns the number of nucleotide difference(s) between two codons.
    '''
    return sum([1 for i in range(len(source)) if source[i] != target[i]])


def find_min_dist(cdn_a, cdn_b):
    dists = [nuc_diff(c1, c2) for c2 in cdn_b for c1 in cdn_a]
    return (min(dists))


def create_min_dist_dictionary(AminoAcid, AA_Codon):
    min_dist_dict = {}

    for i in range(len(AminoAcid)):
        for j in range(len(AminoAcid)):
            a = AminoAcid[i]
            b = AminoAcid[j]
            cdn_a = AA_Codon[a]
            cdn_b = AA_Codon[b]

            min_dist = find_min_dist(cdn_a, cdn_b)
            min_dist_dict[a + b] = min_dist
    return (min_dist_dict)
