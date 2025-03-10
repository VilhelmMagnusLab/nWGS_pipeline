#!/usr/bin/python3
import pandas as pd
from NN_model import NN_classifier

# Retrieve paths from Snakemake's provided input and output variables
model_file_path = snakemake.input.model
bed_file_path = snakemake.input.bed
output_summary_path = snakemake.output.txt
output_votes_path = snakemake.output.votes

# Initialize the model
NN = NN_classifier(model_file_path)
bed_sample = pd.read_csv(bed_file_path, delimiter="\t", index_col=0)  # Load the bed file
#print("NNNNN",bed_sample.shape)
predictions, class_labels, n_features = NN.predict(bed_sample)

# Write predictions to a table
#df = pd.DataFrame({'class': class_labels, 'score': predictions, 'num_features': [n_features]})
df = pd.DataFrame({
    'class': class_labels,
    'score': predictions,
    'num_features': [n_features] * len(class_labels)  # repeat `n_features` if needed
})
df.to_csv(output_votes_path, sep='\t', index=False)
print("XXXXX",df)

# Write summary to a txt file
summary = [
    f'Number of features: {n_features}',
    f'Predicted Class: {class_labels[0]}',
    f'Score: {predictions[0]}'
]

with open(output_summary_path, 'w') as f:
    f.write("\n".join(summary))


# Optional: Uncomment if you need to load and save the BED file
# bed = pd.read_csv(snakemake.input["bed"], delimiter="\t")
# bed.to_csv(snakemake.output['votes'], sep='\t', index=False)



# Read the BED file as per the required format (replace predict_from_bedMethyl)
###bed = pd.read_csv(snakemake.input["bed"], delimiter="\t")
#d.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'coverage', 'MAF']

# Display or use the loaded BED data as needed, but without any prediction or additional processing.
# The following code demonstrates an example of handling the data further, if needed.

# For demonstration purposes, we'll just print the BED file data
###print(bed)

# If there is any further processing or saving needed, you can add it here as appropriate.
# For example:
###bed.to_csv(snakemake.output['votes'], sep='\t', index=False)

# Write summary to txt file
###summary = ['BED file loaded and saved successfully.']

###with open(snakemake.output["txt"], 'w') as f:
   ### f.write("\n".join(summary))
