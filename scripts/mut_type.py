from Bio.SeqUtils import ProtParamData

def analyze_mutation(original, substitute):
    properties = {
        'size': ProtParamData.kd,
        'polarity': ProtParamData.IP,
        'hydrophobicity': ProtParamData.hw
    }

    differences = []
    for prop, values in properties.items():
        # Define a threshold for each property
        threshold = {
            'size': 20,  # Arbitrary example threshold
            'polarity': 1.0,  # Arbitrary example threshold
            'hydrophobicity': 0.5  # Arbitrary example threshold
        }[prop]

        if abs(values[original] - values[substitute]) > threshold:
            differences.append(f"{prop} (original: {values[original]}, substitute: {values[substitute]})")

    if differences:
        return "Non-conservative", differences
    else:
        return "Conservative", []

# Example usage
mutation_type, diff = analyze_mutation('A', 'S')  # Alanine (A) to Serine (S)
print(f"Mutation type: {mutation_type}")
if diff:
    print("Differences:", ", ".join(diff))

