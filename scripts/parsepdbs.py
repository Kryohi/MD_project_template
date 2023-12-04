import sys
import os
import numpy as np
from Bio import PDB, SeqIO, Align

## TODO convert to Biotite which is more complete


def get_pdb_sequence(structure):
    """ Extract sequence of residues from the PDB structure. """
    ppb = PDB.PPBuilder()
    for pp in ppb.build_peptides(structure):
        return pp.get_sequence()
    
def find_gaps(residue_numbers):
    # Find gaps in a sorted list of residue numbers
    gaps = []
    for i in range(1, len(residue_numbers)):
        if residue_numbers[i] != residue_numbers[i-1] + 1:
            gaps.append((residue_numbers[i-1] + 1, residue_numbers[i] - 1))
    return gaps

def determine_chain_type(chain):
    # Simplified logic to determine chain type based on residue names
    for residue in chain.get_residues():
        res_name = residue.get_resname()
        if res_name in PDB.Polypeptide.standard_aa_names:
            return "Protein"
        elif res_name in PDB.Polynucleotide.standard_nt_names:
            return "Nucleic"
    return "Other"


def analyze_secondary_structure(structure):
    """ Analyze the secondary structure composition using DSSP. """
    dssp = PDB.DSSP(structure[0], pdb_path, dssp="dssp")
    ss_counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
    for _, _, ss, _, _, _, _, _, _ in dssp:
        if ss in ss_counts:
            ss_counts[ss] += 1
    return ss_counts

def find_disulfide_bonds(structure):
    disulfide_bonds = []
    cysteines = [residue for residue in structure.get_residues() if residue.get_resname() == 'CYS']

    for i, cys1 in enumerate(cysteines):
        for cys2 in cysteines[i+1:]:
            # Check if the sulfur atoms are within a certain distance (e.g., 2 angstrom is the gromacs default)
            if (cys1['SG'].get_vector() - cys2['SG'].get_vector()).norm() < 2.2:
                disulfide_bonds.append((cys1.get_full_id(), cys2.get_full_id()))

    return disulfide_bonds

def calculate_radius_of_gyration(structure):
    atom_coords = [atom.get_coord() for atom in structure.get_atoms() if atom.get_name() == 'CA']
    if not atom_coords:
        return 'N/A'  # Return 'N/A' if no CA atoms are found

    center_of_mass = np.mean(atom_coords, axis=0)
    rg = np.sqrt(np.mean(np.sum((atom_coords - center_of_mass)**2, axis=1)))

    return rg


def analyze_pdb(pdb_path, fasta_path):
    pdb_id = os.path.basename(pdb_path).split('.')[0]
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    fasta_seq = next(SeqIO.parse(fasta_path, "fasta")).seq

    chain_info = {}
    missing_residues = {}
    ss_composition = {}


    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            chain_type = determine_chain_type(chain)
            chain_info[chain_id] = chain_type

            # Check for gaps in residue numbering
            residue_numbers = [res.get_id()[1] for res in chain if res.get_id()[0] == ' ']
            gaps = find_gaps(residue_numbers)
            if gaps:
                missing_residues[chain_id] = gaps

    # Secondary structure analysis
    #ss_composition = analyze_secondary_structure(structure)

    # Disulfide bonds analysis
    disulfide_bonds = find_disulfide_bonds(structure)

    # Radius of gyration calculation
    radius_of_gyration = calculate_radius_of_gyration(structure)

    # Compile results
    results = {
        "missing loops": missing_residues,
        "secondary_structure": ss_composition,
        "disulfide_bonds": disulfide_bonds,
        "radius_of_gyration": radius_of_gyration
    }
    return results



# Check if at least one argument (PDB file path) is provided
if len(sys.argv) < 2:
    print("Usage: python analyze_pdb.py <pdb_file_path>")
    sys.exit(1)

# Process each PDB file path provided as an argument
pdb_path = sys.argv[1]
fasta_path = pdb_path.replace('.pdb', '.fasta')  # Assuming FASTA file is named similarly
analysis_results = analyze_pdb(pdb_path, fasta_path)
#print(f"Analysis for {pdb_path}:")
print(analysis_results)
print("---------")

# Example usage
#pdb_id = "1IA8"  # Example PDB ID
#pdb_path = f"./pdb_data/{pdb_id}.pdb"  # Replace with your path
#chain_info, missing_residues = analyze_pdb(pdb_path)

