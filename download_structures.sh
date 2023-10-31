#!/usr/bin/env sh

NAME_GENE=$1
PDB_ID=$2
echo

# Base URL for downloading PDB files
PDB_BASE_URL="https://files.rcsb.org/download"

# PDB ID to download.
#PDB_ID="1WMV"

# Download from RCSB PDB if PDB_ID is provided
if [ ! -z "$PDB_ID" ]; then
    if [ ! -f "./00-structures/${PDB_ID}_rcspdb.pdb" ]; then
        curl -o "./00-structures/${PDB_ID}_rcspdb.pdb" "${PDB_BASE_URL}/${PDB_ID}.pdb"
        if [ $? -ne 0 ]; then
            echo "Failed to download ${PDB_ID}_rcspdb.pdb"
            exit 1
        fi
    else
    echo "${PDB_ID}_rcspdb.pdb already exists in 00-structures, skipping download."
    fi
fi


# TODO check foldseek as an alternative source

# Query UniProt for the gene
UNIPROT_RESULT=$(curl -s -X GET "https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(gene:${NAME_GENE})&format=tsv")
UNIPROT_RESULT_JSON=$(curl -s -X GET "https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(gene:${NAME_GENE})&format=json")


# Extract Entry and Entry Name # TODO use json here as well
ENTRY=$(echo "$UNIPROT_RESULT" | awk -F'\t' 'NR==2{print $1}')
ENTRY_NAME=$(echo "$UNIPROT_RESULT" | awk -F'\t' 'NR==2{print $2}')
echo -e "\nFound UniProt entry $ENTRY for the gene $NAME_GENE, Homo Sapiens."
echo "$UNIPROT_RESULT_JSON" > "./00-structures/${ENTRY}_UniProt.json"
jq '.' "./00-structures/${ENTRY}_UniProt.json" > "./00-structures/${ENTRY}_UniProt_pretty.json"
rm "./00-structures/${ENTRY}_UniProt.json"
echo -e "UniProt data saved to ./00-structures/${ENTRY}_UniProt_pretty.json"


# Query AlphaFold using the gene name (homo sapiens version will be requested)
AF_RESULT=$(curl -s -X 'GET' \
  "https://alphafold.com/api/prediction/${ENTRY}" \
  -H 'accept: application/json')

echo "$AF_RESULT" > "./00-structures/${ENTRY}_AF2DB.json"


# Extract the PDB URL from the JSON (always takes the first element)
AF_RESULT_URL=$(echo "$AF_RESULT" | jq -r '.[0].entryId') # could also directly extract pdbUrl etc.
echo -e "Found AlphaFold DB entry $AF_RESULT_URL, now downloading...\n"

if [ ! -z "$AF_RESULT_URL" ]; then
    # Construct the URLs for the AlphaFold structure
    AF_PDB_URL=$(echo "$AF_RESULT" | jq -r '.[0].pdbUrl')
    AF_CIF_URL=$(echo "$AF_RESULT" | jq -r '.[0].cifUrl')
    AF_PAE_URL=$(echo "$AF_RESULT" | jq -r '.[0].paeDocUrl')
    AF_PAEIMG_URL=$(echo "$AF_RESULT" | jq -r '.[0].paeImageUrl')

    # Download the AlphaFold structures
    if [ ! -f "./00-structures/${NAME_GENE}_AF2.pdb" ]; then
        curl -f -o "./00-structures/${NAME_GENE}_AF2.pdb" "$AF_PDB_URL" || { echo "Failed to download ${NAME_GENE}_AF2.pdb"; exit 1; }
        curl -f -o "./00-structures/${NAME_GENE}_AF2.cif" "$AF_CIF_URL" || { echo "Failed to download ${NAME_GENE}_AF2.cif"; exit 1; }
        curl -f -o "./00-structures/${NAME_GENE}_AF2_PAE.json" "$AF_PAE_URL" || { echo "Failed to download ${NAME_GENE}_AF2_PAE.json"; exit 1; }
        curl -f -o "./00-structures/${NAME_GENE}_AF2_PAE.png" "$AF_PAEIMG_URL" || { echo "Failed to download ${NAME_GENE}_AF2_PAE.png"; exit 1; }

        echo "Successfully downloaded data for the ${AF_RESULT_URL} entry!"
    else
        echo "${NAME_GENE}_AF2.pdb already exists in 00-structures, skipping download."
    fi
else
    echo "No AlphaFold result found for UniProt entry ${ENTRY}"
fi

echo

conda env config vars set GENE_NAME=$NAME_GENE
conda env config vars set UNIPROT_ID=$ENTRY
conda env config vars set PDB_FILE="${NAME_GENE}_AF2.pdb"
echo "Environmental variable UNIPROT_ID set for the environment $CONDA_DEFAULT_ENV"


echo
