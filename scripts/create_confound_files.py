import os
import json
import pandas as pd
import numpy as np
import logging

def read_tsv(inTSV):
	dataTSV = pd.read_table(inTSV)
	dataTSV.head()
	return dataTSV

def read_json(inJSON):
	with open(inJSON,'rt') as cj:
		dataJSON = json.load(cj)
	return dataJSON

logging.basicConfig(filename=snakemake.log.logfile,filemode='w',format='%(name)s - %(levelname)s - %(message)s',level=logging.INFO)

logging.info('Test logging message')

# Create empty confound_file for case of no experimental confounds
if (snakemake.params.confound_name == 'noconfounds'):
	open(snakemake.output.confound_file,'a').close()
else:
	# Read in confounds.tsv from fmriprep and confounds dictionary from config.yaml
	confounds = read_tsv(snakemake.input.confounds_tsv)
	confounds_dict = read_json(snakemake.input.confounds_dictionary)

	# Extract confound set based on confound_name parameter
	# Save as text file for input into FSL design.fsf
	confound_set = confounds_dict.get(snakemake.params.confound_name)
	confound_set_save = confounds[confound_set]
	confound_set_save = np.nan_to_num(np.array(confound_set_save))
	np.savetxt(snakemake.output.confound_file, confound_set_save, delimiter='\t')
