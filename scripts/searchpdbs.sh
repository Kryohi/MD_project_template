#!/bin/bash

# Define the UniProt API endpoints for searching and fetching entry data
UNIPROT_SEARCH_ENDPOINT="https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_ENTRY_ENDPOINT="https://rest.uniprot.org/uniprotkb"

# Define the gene name and organism ID
gene_name=$1
organism_id=$2 

# Create a directory to save the data
data_dir="$PROJECT_DIR/data/00-structures"
#mkdir -p "$data_dir"

# Fetch entry data for the specified gene and organism from UniProt
# Use TSV format for the request
query="gene_exact:${gene_name} AND organism_id:${organism_id}"
tsv_data=$(curl -s -G --data-urlencode "query=${query}" --data-urlencode "format=tsv" "$UNIPROT_SEARCH_ENDPOINT")

# Check if any data was found in TSV
if [ -z "$tsv_data" ]; then
    echo "No data found for the gene: $gene_name and organism ID: $organism_id"
    exit 1
fi

# Extract the first UniProt Entry from the TSV data
entry_name=$(echo "$tsv_data" | awk 'NR==2 {print $1}')

# Fetch the JSON data for the extracted entry name from UniProt
json_file_path="${data_dir}/${entry_name}.json"
curl -s "${UNIPROT_ENTRY_ENDPOINT}/${entry_name}.json" -o "$json_file_path"

# Check if the JSON file was downloaded successfully
if [ ! -f "$json_file_path" ]; then
    echo "Failed to download JSON data for the entry name: $entry_name"
    exit 1
fi

jq '.' "$json_file_path" > "$json_file_path.1" # prettify
rm "$json_file_path"
mv "$json_file_path.1" "$json_file_path"

echo
echo "JSON data saved to: $json_file_path"
echo

# Extract PDB IDs from the JSON file
pdb_ids=$(jq -r '.uniProtKBCrossReferences[] | select(.database == "PDB") | .id' "$json_file_path" | sort | uniq)

# Check if any PDB IDs were extracted
if [ -z "$pdb_ids" ]; then
    echo "No PDB IDs found in the JSON file."
    exit 1
fi

echo "Extracted PDB IDs:"
# Print the extracted PDB IDs in a single line, separated by commas
echo "$pdb_ids" | tr '\n' ', ' | sed 's/,$//'
echo
echo

# RCSB PDB API base URL
RCSB_API_BASE="https://data.rcsb.org/rest/v1/core/entry"


# Initialize counters
total_pdbs=$(echo "$pdb_ids" | wc -l)
downloaded_pdbs=0

echo "Starting to download PDB JSON files..."

# Loop through PDB IDs and fetch data from RCSB PDB
for pdb_id in $pdb_ids; do
    # Fetch data from RCSB PDB for the given PDB ID
    pdb_file_path="${data_dir}/${pdb_id}.json"
    if (curl -s "${RCSB_API_BASE}/${pdb_id}" -o "$pdb_file_path"); then
        if [ ! -s "$pdb_file_path" ]; then
            echo "Warning: Empty JSON file for PDB ID $pdb_id"
            continue  # Skip processing this file
        fi
        jq '.' "$pdb_file_path" > "$pdb_file_path.1" # prettify
        rm "$pdb_file_path"
        mv "$pdb_file_path.1" "$pdb_file_path"
        ((downloaded_pdbs++))
    else
        echo "Failed to download data for PDB ID $pdb_id"
    fi
done

echo "Downloaded $downloaded_pdbs out of $total_pdbs PDB JSON files."


# File to save the table in CSV format
table_file="${data_dir}/pdb_table.csv"

# Print table header to the file
printf "Index,PDB ID,Method,Resolution,Modeled Residues,Deposited Residues,Has Nucleic,Has Ligands,Title\n" > "$table_file"

# Initialize row counter
row_index=0

# Loop through the downloaded JSON files
for pdb_file in "${data_dir}"/*.json; do
    ((row_index++))  # Increment row counter
    pdb_id=$(basename "$pdb_file" .json)

    # Extract properties from the JSON file
    method=$(jq -r '.exptl[0].method // "N/A"' "$pdb_file")
    resolution=$(jq -r '.reflns[0].d_resolution_high // "N/A"' "$pdb_file")
    modeled_residues=$(jq -r '.rcsb_entry_info.deposited_modeled_polymer_monomer_count // "N/A"' "$pdb_file")
    deposited_residues=$(jq -r '.rcsb_entry_info.deposited_polymer_monomer_count // "N/A"' "$pdb_file")
    has_nucleic=$(jq -r '[.refine_hist[] | .pdbx_number_atoms_nucleic_acid | select(. != null) | tonumber] | any > 0' "$pdb_file")
    has_ligands=$(jq -r '[.refine_hist[] | .pdbx_number_atoms_ligand | select(. != null) | tonumber] | any > 0' "$pdb_file")
    title=$(jq -r '.struct.title // "N/A"' "$pdb_file")

    # Convert true/false to 1/0 for has_nucleic and has_ligands
    has_nucleic=$( [[ $has_nucleic == true ]] && echo "1" || echo "0" )
    has_ligands=$( [[ $has_ligands == true ]] && echo "1" || echo "0" )

    # Add a line to the table file
    printf "%d,%s,%s,%s,%s,%s,%s,%s,%s\n" "$row_index" "$pdb_id" "$method" "$resolution" "$modeled_residues" "$deposited_residues" "$has_nucleic" "$has_ligands" "$title" >> "$table_file"
done

echo "Table of PDB properties saved to: $table_file"

# Display the CSV file as a table
echo "Displaying the CSV file as a table:"
column -t -s ',' "$table_file" | less -S

echo "End of table."


# Ask the user for a list of indices
read -p "Enter a list of indices separated by space that you want to download: " -a indices

# Path to the CSV file
csv_file="${data_dir}/pdb_table.csv"

# RCSB PDB file download endpoints
PDB_DOWNLOAD_BASE="https://files.rcsb.org/download"
FASTA_DOWNLOAD_BASE="https://www.rcsb.org/fasta/entry"

# Loop over the provided indices to fetch the corresponding PDB IDs
for index in "${indices[@]}"; do
    echo
    # Extract the PDB ID corresponding to the index
    # Note: CSV is 1-indexed, and awk is 1-indexed
    pdb_id=$(awk -v idx="$index" -F ',' 'NR==idx+1 {print $2}' "$csv_file")

    # Check if pdb_id was found
    if [ -z "$pdb_id" ]; then
        echo "No PDB ID found for index: $index"
        continue
    fi

    # Download PDB file
    pdb_file="${data_dir}/${pdb_id}.pdb"
    curl -s "${PDB_DOWNLOAD_BASE}/${pdb_id}.pdb" -o "$pdb_file"

    # Download mmCIF file
    mmcif_file="${data_dir}/${pdb_id}.cif"
    curl -s "${PDB_DOWNLOAD_BASE}/${pdb_id}.cif" -o "$mmcif_file"

    # Download FASTA file
    fasta_file="${data_dir}/${pdb_id}.fasta"
    curl -s "${FASTA_DOWNLOAD_BASE}/${pdb_id}" -o "$fasta_file"

    echo "Downloaded PDB, mmCIF and fasta files for PDB ID: $pdb_id in $data_dir"


    if [ -f "$pdb_file_path" ]; then
        python parsepdbs.py "$pdb_file"
    else
        echo "PDB file not found for ID: $pdb_id"
    fi
    echo

    # TODO Try to download the paper of the PDB
    # Path to the already downloaded JSON file for metadata
    metadata_file="${data_dir}/${pdb_id}.json"

    # Check if the JSON metadata file exists
    if [ -f "$metadata_file" ]; then
        # Extract DOI from the metadata JSON file
        doi=$(jq -r '.rcsb_primary_citation.pdbx_database_id_doi // "N/A"' "$metadata_file")

        if [ "$doi" != "N/A" ]; then
            echo "DOI for PDB ID $pdb_id: $doi"
            # Construct the URL to the paper (Note: This might not lead directly to a PDF)
            paper_url="https://doi.org/$doi"
            echo "URL for paper: $paper_url"
        else
            echo "No DOI found for PDB ID $pdb_id"
        fi
    else
        echo "Metadata JSON file not found for PDB ID $pdb_id"
    fi

done
echo

# Delete all the PDB JSON files but keep the UniProt entry JSON file
rm "${data_dir}"/[0-9]*.json


