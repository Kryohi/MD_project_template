 
import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter

# Load the source PDB
u_source = mda.Universe('../1yfh_1.pdb')

# Load the destination PDB
u_dest = mda.Universe('../AF-P16455-F1-aligned.pdb')

# Select atoms.
selected_atoms = u_source.select_atoms('resname XCY or name ZN or nucleic')


# Count the number of selected atoms
num_atoms = len(selected_atoms)
print(f"Number of selected atoms: {num_atoms}")

# Print the details of residues and chains
residues = set()  # Using a set to avoid duplicates
chains = set()  # Using a set to avoid duplicates

for atom in selected_atoms:
    residues.add(atom.residue)
    chains.add(atom.segment)

# Print the residues
print("Residues containing selected atoms:")
#for residue in residues:
    #print(f"Residue: {residue.name}, Residue ID: {residue.resid}")

# Print the chains
print("Chains containing selected atoms:")
for chain in chains:
    print(f"Chain ID: {chain.segid}")


# Combine the destination atoms and the selected atoms into a new Universe
all_atoms = mda.Merge(u_dest.atoms, selected_atoms)

# Write the new Universe to a PDB file
with PDBWriter('../AF2_ZN_nucleic_XCY.pdb') as writer:
    writer.write(all_atoms)