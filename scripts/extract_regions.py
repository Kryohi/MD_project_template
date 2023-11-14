import json
import os


def extract_regions_from_uniprot():

    # Retrieve the UNIPROT_ID environment variable
    uniprot_id = os.environ.get('UNIPROT_ID')

    # Check if the UNIPROT_ID is set
    if uniprot_id is not None:
        os.getcwd()
        # Construct the file path
        file_path = f"../data/00-structures/external/{uniprot_id}_UniProt_pretty.json"
        print("File path:", file_path)
    else:
        print("UNIPROT_ID is not set")
        exit(1)


    # Load the JSON data
    with open(file_path, 'r') as file:
        uniprot_data = json.load(file)


    # Extracting relevant information
    protein_description = uniprot_data['results'][0].get('proteinDescription', {})
    genes = uniprot_data['results'][0].get('genes', [])
    comments = uniprot_data['results'][0].get('comments', [])
    features = uniprot_data['results'][0].get('features', [])
    sequence = uniprot_data['results'][0].get('sequence', {}).get('sequence', '')


    # Handling cases where 'type' key is missing and providing a default value
    comments_by_type = {}
    for comment in comments:
        comment_type = comment.get('type', 'Unknown')
        if comment_type not in comments_by_type:
            comments_by_type[comment_type] = []
        comments_by_type[comment_type].append(comment)

    features_by_type = {}
    for feature in features:
        feature_type = feature.get('type', 'Unknown')
        if feature_type not in features_by_type:
            features_by_type[feature_type] = []
        features_by_type[feature_type].append(feature)

    # Displaying a summary of each section to understand their content
    summary = {
        "proteinDescription": protein_description,
        "genes": genes,
        "comments": comments_by_type,
        "features": features_by_type,
        "sequence_length": len(sequence),
        "sequence_preview": sequence[:50] + "..."  # Displaying the first 50 characters of the sequence
    }

    #summary


    # Extracting domain and region information along with their residue intervals
    domains = []
    regions = []

    for feature in features:
        feature_type = feature.get('type')
        if feature_type == 'Domain':
            domains.append({
                'name': feature.get('description'),
                'start': feature['location']['start']['value'],
                'end': feature['location']['end']['value']
            })
        elif feature_type == 'Region':
            regions.append({
                'name': feature.get('description'),
                'start': feature['location']['start']['value'],
                'end': feature['location']['end']['value']
            })


    # Return domains and regions
    return domains, regions


# Allow the script to be run standalone or imported as a module
if __name__ == "__main__":
    domains, regions = extract_regions_from_uniprot()
    print("Domains:", domains)
    print("Regions:", regions)

