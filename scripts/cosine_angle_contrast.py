import os
import json
from nilearn.image import load_img 
import numpy as np

contrast_vector = snakemake.params.contrast_vector

for file in range(len(sublist)):
	if num_evs == 1:
		peFile = load_img(sublist[file])
		peData = peFile.get_data()
		if peData.ndim == 3:
			peData = peData
		elif peData.ndim == 4:
			peData = peData[:,:,:,0]
		peData2Dtemp = np.reshape(peData,((np.shape(peData)[0]*np.shape(peData)[1]*np.shape(peData)[2]),-1))
		peData2D.append(peData2Dtemp)
	else:
		peFile = load_img(sublist[file])
		peData = peFile.get_data()
		peData = peData[:,:,:,0:numEVs]
		peData2Dtemp = np.reshape(peData,((np.shape(peData)[0]*np.shape(peData)[1]*np.shape(peData)[2]),np.shape(peData)[3]))
		peData2D.append(peData2Dtemp)
			
	if numEVs == 1:
		peData2D = np.squeeze(np.array(peData2D).T)
	else:
		peData2D = np.array(peData2D).T
		peData2D = np.swapaxes(peData2D,0,1)
			
	peDataStack = peData2D
	

	for run in range(numRuns):
		testing = peDataStack[:,:,run]
		training = np.mean(np.c_[peDataStack[:,:,0:run],peDataStack[:,:,run+1:numRuns]],axis=2)
		numerator = np.dot(contrast_vector.T,np.dot(training.T,np.dot(testing,contrast_vector)))
		denominator = np.sqrt(np.dot(contrast_vector.T,np.dot(training.T,np.dot(training,contrast_vector)))*np.dot(contrast_vector.T,np.dot(testing.T,np.dot(testing,contrast_vector))))
		temp = numerator/denominator
		cosine_temp.append(np.squeeze(temp))

    mean_cosine.append(np.mean(np.array(cosine_temp)))
	cosine.append(cosine_temp)
			

    cosine_dictionary = dict(zip(snakemake.params.confound_names,cosine))
	with open(snakemake.output.cosine_dictionary, 'w') as fp:
		json.dump(cosine_dictionary, fp, indent=4, separators=(',', ': '))

	mean_cosine_dictionary = dict(zip(snakemake.params.confound_names,mean_cosine))
	with open(snakemake.output.mean_cosine_dictionary, 'w') as fp:
		json.dump(mean_cosine_dictionary, fp, indent=4, separators=(',', ': '))