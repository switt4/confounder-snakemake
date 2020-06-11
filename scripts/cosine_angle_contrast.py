import os
import json
from nilearn.image import load_img 
import numpy as np

sublist = snakemake.input.stats_files
num_runs = len(sublist)

num_evs = len(snakemake.params.trial_names)
contrast_vector = snakemake.params.contrast_vector

contrast_vector = np.fromstring(contrast_vector,dtype=float,sep=',')
contrast_vector = contrast_vector[...,np.newaxis]

confound_names = snakemake.params.confound_name

mean_cosine = []
cosine_temp = []
cosine =[]

# by keeping the EVs dim of in dim 2 as it is with multiple EVs this section can be significatly shortened
for file in range(len(sublist)):
	peFile = load_img(sublist[file])
	peData = peFile.get_data()
	
	if peData.ndim == 4:
		peData = peData[:,:,:,0:num_evs]
	peData2Dtemp = np.reshape(peData,(np.shape(peData)[0]*np.shape(peData)[1]*np.shape(peData)[2],-1))

	#this is what allows for the dim 2 of one, intialize with new axis then concatenate along it
	if file>0:
		peData2D = np.concatenate([peData2D, peData2Dtemp[:,:,np.newaxis]], axis=2)
	else:
		peData2D = peData2Dtemp[:,:,np.newaxis]

for run in range(num_runs):
	testing = peData2D[:,:,run]
	training = np.mean(np.c_[peData2D[:,:,0:run],peData2D[:,:,run+1:num_runs]],axis=2)
	# train*test/((train*train)*(test*test))
	# Trace(a dot B)==Trace(a*B) was used to simplify these equations
	numerator = np.dot(contrast_vector.T,np.dot(training.T,np.dot(testing,contrast_vector)))
	denominator = np.sqrt(np.trace(contrast_vector.T*np.dot(training.T,np.dot(training,contrast_vector)))*np.trace(contrast_vector.T*np.dot(testing.T,np.dot(testing,contrast_vector))))
	temp = numerator/denominator
	temp = np.squeeze(temp)
	temp = temp.tolist()
	cosine_temp.append(temp)

mean_cosine.append(np.mean(np.array(cosine_temp)))
cosine.append(cosine_temp)

# Create dictionary of cosine similarity values and confound_name
cosine_dictionary = dict(zip(confound_names[:],cosine))
cosine_dictionary[confound_names[:]] = cosine_dictionary.pop(confound_names[0])
with open(snakemake.output.cosine_dictionary, 'w') as fp:
	json.dump(cosine_dictionary, fp, indent=4, separators=(',', ': '))

mean_cosine_dictionary = dict(zip(confound_names[:],mean_cosine))
mean_cosine_dictionary[confound_names[:]] = mean_cosine_dictionary.pop(confound_names[0])
with open(snakemake.output.mean_cosine_dictionary, 'w') as fp:
	json.dump(mean_cosine_dictionary, fp, indent=4, separators=(',', ': '))