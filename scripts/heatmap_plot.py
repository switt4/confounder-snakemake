import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

def read_tsv(inTSV):
	dataTSV = pd.read_table(inTSV)
	dataTSV.head()
	return dataTSV

def read_json(inJSON):
	with open(inJSON,'rt') as cj:
		dataJSON = json.load(cj)
	return dataJSON

def heatmap_plot(inArray,XLabels,YLabels,Title,Filename):
	fig,ax = plt.subplots()
	im = ax.imshow(inArray)
	ax.set_title(Title,fontsize='xx-large')
	ax.set_xticks(np.arange(np.shape(inArray)[1]))
	ax.set_yticks(np.arange(np.shape(inArray)[0]))
	ax.set_xticklabels(XLabels,rotation=45,ha='right',rotation_mode='anchor')
	ax.set_yticklabels(YLabels)
	cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.35)
	cbar.ax.set_xlabel('correlation')
	cbar.set_clim(vmin=-0.3, vmax=0.3)
	plt.tight_layout()
	plt.savefig(Filename,format='svg')
	plt.close()

confounds = read_tsv(snakemake.input.confounds_tsv)

confounds_dictionary = read_json(snakemake.input.confounds_dictionary)
confounds_names = list(confounds_dictionary.values())
flat_names = [item for sublist in confounds_names for item in sublist]
flat_names_unique = list(set(flat_names))
flat_names_unique.sort()

bold_signal = np.loadtxt(snakemake.input.bold_signal_file)

bold_confound_correlation = []

for bold in range(len(snakemake.params.trial_names)):
    for conf in range(len(flat_names_unique)):
        corr_temp = np.corrcoef(bold_signal[:,bold],np.nan_to_num(np.array(confounds[flat_names_unique[conf]])))
        bold_confound_correlation.append(corr_temp[0,1])

bold_confound_correlation = np.array(bold_confound_correlation)
bold_confound_correlation = np.reshape(bold_confound_correlation,(len(snakemake.params.trial_names),len(flat_names_unique)))

heatmap_plot(bold_confound_correlation,flat_names_unique,snakemake.params.trial_names,'Correlation between Bold Signal and Confounds',snakemake.output[0])