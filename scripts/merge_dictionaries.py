import json

cosine_dictionaries = snakemake.input.cosine_dictionaries

# Concatenate individual dictionary files into a list of dictionaries
merged_temp = []
for f in cosine_dictionaries:
    with open(f, "r") as infile:
        merged_temp.append(json.load(infile))

# Convert list of dictionaries into a single dictionary
merged_dictionary = dict((key,d[key]) for d in merged_temp for key in d)

with open(snakemake.output.merged_cosine_dictionary, "w") as outfile:
     json.dump(merged_dictionary, outfile, indent=4, separators=(',', ': '))