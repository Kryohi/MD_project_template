#!/usr/bin/env sh

# Define the gene of interest
gene_name="WWOX"

cd ./data/00-metadata


echo "Gene ID:"
curl -X POST -H "Content-Type: application/xml" -d "<?xml version="1.0" encoding="UTF-8"?>
<\!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "external_gene_name" value = $gene_name/>
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "description" />
        <Attribute name = "phenotype_description" />
    </Dataset>
</Query>' 'https://www.ensembl.org/biomart/martservice/results" #> metadata.tsv



# Query BioMart for the canonical transcript of the specified gene
curl -X POST -H "Content-Type: application/xml" -d "<?xml version='1.0' encoding='UTF-8'?>
<\!DOCTYPE Query>
<Query virtualSchemaName='default' formatter='TSV' header='0' uniqueRows='0' count='' datasetConfigVersion='0.6'>
    <Database name='ensembl_mart_110'>
        <Dataset name='hsapiens_gene_ensembl' interface='default'>
            <Attribute name='ensembl_gene_id' />
            <Attribute name='ensembl_transcript_id' />
            <Attribute name='transcript_biotype' />
            <Attribute name='canonical_transcript' />
            <Attribute name='protein_id' />
            <Attribute name='external_gene_name' />
            <Filter name='external_gene_name' value='${gene_name}'/>
        </Dataset>
    </Database>
</Query>" "https://www.ensembl.org/biomart/martservice/results" > wwox_canonical_transcript.tsv

awk '$4 == "1"' wwox_canonical_transcript.tsv > $gene_name_only_canonical.tsv
