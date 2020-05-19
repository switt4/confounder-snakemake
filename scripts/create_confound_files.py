import os
import json
import pandas as pd
import numpy as np

def read_tsv(inTSV):
	dataTSV = pd.read_table(inTSV)
	dataTSV.head()
	return dataTSV

def read_json(inJSON):
	with open(inJSON,'rt') as cj:
		dataJSON = json.load(cj)
	return dataJSON

# Read in confounds.tsv from fmriprep and confounds dictionary from config.yaml
confounds = read_tsv(snakemake.input.confounds_tsv)
confounds_dict = read_json(snakemake.input.confounds_dictionary)

# check if confound name is 'none', skip making confound file if true
if (snakemake.params.confound_name == 'none'):
    continue
else:
    # Extract confound set based on confound_name parameter
    # Save as text file for input into FSL design.fsf
    confound_set = confounds_dict.get(snakemake.params.confound_name)
    confound_set_save = confounds[confound_set]
    confounds_set_save = np.nan_to_num(np.array(confounds_set_save))
    np.savetxt(snakemake.output.confound_file, confounds_set_save, delimiter='\t')

