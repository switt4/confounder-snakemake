#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:47:03 2019

@author: Suzanne T. Witt, Ph.D.
"""

import os
import sys
import re
import fnmatch
import json
import pandas as pd
import argparse
import subprocess
import nilearn
import numpy as np
import xml.etree.ElementTree as ET
from glob import glob
from nilearn.image import (load_img, resample_to_img)
import fsl
import fsl.data.featanalysis as FA
import matplotlib.pyplot as plt

def run(command, env={}):
	merged_env = os.environ
	merged_env.update(env)
	process = subprocess.Popen(command, stdout=subprocess.PIPE,
							   stderr=subprocess.STDOUT, shell=True,
							   env=merged_env)
	while True:
		line = process.stdout.readline()
		line = str(line, 'utf-8')[:-1]
		print(line)
		if line == '' and process.poll() != None:
			break
	if process.returncode != 0:
		raise Exception("Non zero return code: %d"%process.returncode)

def read_tsv(inTSV):
	dataTSV = pd.read_table(inTSV)
	dataTSV.head()
	#confoundsHeaders = list(confounds)
	return dataTSV

def read_json(inJSON):
	with open(inJSON,'rt') as cj:
		dataJSON = json.load(cj)
	return dataJSON

def parse_compcor_json(inDict,mask):
	maskCompCor = []
	for key, subdict in inDict.items():
		sublist = list(subdict.values())
		if mask in sublist:
			maskCompCor.append(key)
	return maskCompCor

def confound_plot(timeScale,inConfound,filename,yScale=[],legend=[]):
	fig,ax = plt.subplots()
	if len(np.shape(inConfound)) == 1:
		ax.plot(timeScale,inConfound)
	else:
		for c in range(np.shape(inConfound)[1]):
			ax.plot(timeScale,inConfound[:,c])
	ax.legend(legend,loc='upper left')
	ax.set_xlabel('time (sec)')
	ax.set_ylabel(yScale)
	plt.tight_layout()
	plt.savefig(filename,format='svg')
	plt.close()


def barh_plot(inDict,Title,filename):
	fig,ax = plt.subplots()
	temp = sorted((value,key) for (key,value) in inDict.items())
	sortedDict = dict([(key, value) for (value, key) in temp])
	ax.barh(list(sortedDict.keys()),list(sortedDict.values()),align='center')
	ax.set_title(Title)
	plt.tight_layout()
	plt.savefig(filename,format='svg')
	plt.close()

parser = argparse.ArgumentParser(description='Example BIDS App entrypoint script.')
parser.add_argument('bids_dir', help='The BIDS directory of the input dataset '
					'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
					'should be stored. If you are running group level analysis '
					'this folder should be prepopulated with the results of the'
					'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
					'Multiple participant level analyses can be run independently '
					'(in parallel) using the same output_dir.',
					choices=['participant'])
parser.add_argument('--fmriprep_dir', help='Specify the same output directory as when '
					'fmriprep was run on the dataset.  E.g., there should be a sub-directory '
					'called, "fmriprep" that was created when fmriprep ran.'
					'The script assumes that the fmriprep output directory is '
					'still formatted according to the standard fmriprep output directory.')
parser.add_argument('--feat_dir', help='Specify the directory where the sample design.fsf '
					'file is saved.')
parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
				   'corresponds to sub-<participant_label> from the BIDS spec '
				   '(so it does not include "sub-"). If this parameter is not '
				   'provided all subjects should be analyzed. Multiple '
				   'participants can be specified with a space separated list.',
				   nargs = "+")
parser.add_argument('--task_label', help='Enter the "task-<task_label>" label for the '
					'task you wish to analyze.  All runs of this task '
					'must have a valid "events.tsv" files.',
					nargs = "+")
parser.add_argument('--HarvardOxford_region', help='Specify region from Harvard Oxford cortical atlas to test. '
					'Use the value in $FSLDIR/data/atlases/HarvardOxford-Cortical.xml; '
					'correction for values starting at 0 are applied automatically. '
					'If parameter is not set, the "Intracalcarine_cortex" (label = 23) will '
					'be tested as default.',
					nargs="+")
#parser.add_argument('-v', '--version', action='version',
#					version='BIDS-App example version {}'.format(__version__))


args = parser.parse_args()

if args.fmriprep_dir is None:
	sys.exit('Error: You must specify an fmnriprep directory.')

if args.feat_dir is None:
	sys.exit('Error: You must specify a directory containing a sample design.fsf file.')

if args.task_label is None:
	sys.exit('Error: You must specify a single task to analyze.')
else:
	taskLabel = str(args.task_label[0])

if not os.path.exists(os.path.join(args.feat_dir,"task-%s"%taskLabel)):
	os.makedirs(os.path.join(args.feat_dir,"task-%s"%taskLabel))

# load in Harvard Oxford atlases and prepare for region extraction
atlasCortexLabels = []
atlasCortexFile = os.path.join('$FSLDIR','data','atlases','HarvardOxford-Cortical.xml')
atlasCortex = ET.parse(os.path.expandvars(atlasCortexFile))
atlasCortexRoot = atlasCortex.getroot()
for roi in range(len(atlasCortexRoot[1])):
	atlasCortexLabels.append(atlasCortexRoot[1][roi].text)

# get full path to $FSLDIR
FSLDIR = os.path.expandvars('$FSLDIR')

# get selected region from Harvard Oxford atlas
if args.HarvardOxford_region:
	testRegionNumber = atlasCortexLabels.index(args.HarvardOxford_region) + 1
else:
	testRegionNumber = atlasCortexLabels.index('Intracalcarine Cortex') + 1

roiRegionNumber = testRegionNumber - 1

subjectsToAnalyze = []
# only for a subset of subjects
if args.participant_label:
	subjectsToAnalyze = args.participant_label
# for all subjects
else:
	subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
	subjectsToAnalyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

if not os.path.exists(os.path.join(args.output_dir,"temp")):
	os.makedirs(os.path.join(args.output_dir,"temp"))

for subjectLabel in subjectsToAnalyze:
	if not os.path.exists(os.path.join(args.output_dir,"sub-%s"%subjectLabel)):
		os.makedirs(os.path.join(args.output_dir,"sub-%s"%subjectLabel))
	if not os.path.exists(os.path.join(args.output_dir,"sub-%s"%subjectLabel,"temp")):
		os.makedirs(os.path.join(args.output_dir,"sub-%s"%subjectLabel,"temp"))
	if not os.path.exists(os.path.join(args.output_dir,"sub-%s"%subjectLabel,"figures")):
		os.makedirs(os.path.join(args.output_dir,"sub-%s"%subjectLabel,"figures"))
	if not os.path.exists(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel)):
		os.makedirs(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel))

# create roi
outROIFilename = "tpl-MNI152NLin2009cAsym_res-02_atlas-HOCPA_desc-th0_dseg-%d.nii.gz"%(roiRegionNumber)
outFileROI = os.path.join(args.output_dir,"temp",outROIFilename)
cmdROI = "fslmaths /Users/switt/Downloads/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_atlas-HOCPA_desc-th0_dseg.nii.gz -thr %d -uthr %d -bin %s"%(testRegionNumber,testRegionNumber,outFileROI)
print(cmdROI)
run(cmdROI)

# define empty lists to save RSS and TSS to until I figure out something better
rssRun = []
tssRun = []

# running participant level
if args.analysis_level == "participant":

	# find all func files and calculate effect of denoising them
	for subjectLabel in subjectsToAnalyze:
		print('starting subject:%s'%subjectLabel)
		funcFileList = glob(os.path.join(args.fmriprep_dir,"fmriprep","sub-%s"%subjectLabel,"func","*task-%s*_bold.nii*"%args.task_label[0])) + glob(os.path.join(args.fmriprep_dir,"fmriprep","sub-%s"%subjectLabel,"ses-*","func","*task-%s*_bold.nii*"%args.task_label[0]))
		for funcFile in funcFileList:
			print('working on functional run: %s'%os.path.split(funcFile)[-1])
			residualsVariance = []
			dataVariance = []
			stringMatch = re.search(r'run+-\d+',funcFile)
			runString = stringMatch.group(0)
			confoundTSVFilename = os.path.split(funcFile)[-1].replace("space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz","desc-confounds_regressors.tsv")
			confoundTSV = os.path.join(os.path.split(funcFile)[0],confoundTSVFilename)
			confounds = read_tsv(confoundTSV)

			bidsDir = os.path.split(funcFile)[0].replace(os.path.join(args.fmriprep_dir,"fmriprep"),args.bids_dir)
			funcJSONFilename = os.path.split(funcFile)[-1].replace("space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz","bold.json")
			funcJSONFile = os.path.join(bidsDir,funcJSONFilename)
			with open(funcJSONFile,'rt') as fj:
				taskInfo = json.load(fj)
			TR = taskInfo['RepetitionTime']
			print('TR = %d sec'%TR)

			funcImg = load_img(funcFile)
			numBOLDScans = funcImg.header['dim'][4]
			print('Number of BOLD scans: %d'%numBOLDScans)
			timeScale = np.arange(numBOLDScans) 

			roiResampleFilename = os.path.split(outFileROI)[-1].replace("desc-th0_","desc-th0-resample_")
			roiResampleFile = os.path.join(args.output_dir,"temp",roiResampleFilename)
			roiImg = load_img(outFileROI)
			print('Resampling %s to %s'%(os.path.split(outFileROI)[-1],os.path.split(funcFile)[-1]))
			resampledROIImg = resample_to_img(roiImg,funcImg)
			print('Saving resampled file: %s'%roiResampleFilename)
			resampledROIImg.to_filename(roiResampleFile)

			# load in design.fsf
			evNames = []
			design = FA.loadSettings(args.feat_dir)
			for ev in range(int(design['evs_orig'])):
				evTemp = ev + 1
				evNames.append(design['evtitle%d'%evTemp])

			# load in events.tsv file
			eventsFilename = os.path.split(funcFile)[-1].replace("space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz","events.tsv")
			print('Events file: %s'%eventsFilename)
			eventsFile = os.path.join(bidsDir,eventsFilename)
			taskTimings = read_tsv(eventsFile)

			for ev in range(len(evNames)):
				evTemp = ev + 1
				timings = taskTimings[taskTimings.trial_type == evNames[ev]]
				timings = timings[['onset','duration']]
				timings = np.column_stack((timings,np.ones(len(timings))))
				timingsFilename = eventsFilename.replace("events.tsv","onsets-%s.txt"%(evNames[ev]))
				timingsFile = os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,timingsFilename)
				np.savetxt(timingsFile,timings,delimiter='\t')  
				custom = [('custom%d'%evTemp,'"%s"'%timingsFile)]
				design.update(custom)

			# first step in tedious process of rebuilding design.fsf 
			feat_files = [('feat_files(1)','"%s"'%funcFile)]
			design.update(feat_files)
			npts = [('npts','%d'%numBOLDScans)]
			design.update(npts)

			featwatcher = [('featwatcher_yn',0)]
			design.update(featwatcher)
			regstandard = [('regstandard','"%s/data/standard/MNI152_T1_2mm_brain"'%FSLDIR)]
			design.update(regstandard)
			confoundevs = [('confoundevs',1)]
			design.update(confoundevs)
			gdc = [('gdc', '""')]
			design.update(gdc)
			motionevsbeta = [('motionevsbeta','""')]
			design.update(motionevsbeta)
			scriptevsbeta = [('scriptevsbeta','""')]
			design.update(scriptevsbeta)
			alternative_mask = [('alternative_mask','""')]
			design.update(alternative_mask)
			init_initial_highres = [('init_initial_highres','""')]
			design.update(init_initial_highres)
			init_highres = [('init_highres','""')]
			design.update(init_highres)
			init_standard = [('init_standard','""')]
			design.update(init_standard)

			# clunky way to create confound matrices for fsl_glm
			confoundsToTest = []
			confoundsToTestNames = []
			print('Extracting confounds from: %s'%confoundTSVFilename)
			if os.path.isfile(confoundTSV.replace(".tsv",".json")):
				confoundJSONFilename = os.path.split(confoundTSV)[-1].replace(".tsv",".json")
				confoundJSON = os.path.join(os.path.split(funcFile)[0],confoundJSONFilename)
				compcorJSON = read_json(confoundJSON)

				combinedMask = parse_compcor_json(compcorJSON,"combined")
				sixCombinedACompCor = confounds[combinedMask[0:5]]

				wmMask = parse_compcor_json(compcorJSON,"WM")
				sixWMACompCor = confounds[wmMask[0:5]]

				csfMask = parse_compcor_json(compcorJSON,"CSF")
				sixCSFACompCor = confounds[csfMask[0:5]]

				sixTCompCor = confounds[['t_comp_cor_00','t_comp_cor_01','t_comp_cor_02','t_comp_cor_03','t_comp_cor_04','t_comp_cor_05']]

				sixMotion = confounds[['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']]
				twelveMotion = confounds[['trans_x','trans_x_derivative1','trans_y','trans_y_derivative1','trans_z','trans_z_derivative1',
										 'rot_x','rot_x_derivative1','rot_y','rot_y_derivative1','rot_z','rot_z_derivative1']]
				twentyfourMotion = confounds[['trans_x', 'trans_x_derivative1', 'trans_x_derivative1_power2', 'trans_x_power2', 
											'trans_y', 'trans_y_derivative1', 'trans_y_derivative1_power2', 'trans_y_power2', 
											'trans_z', 'trans_z_derivative1', 'trans_z_power2', 'trans_z_derivative1_power2', 
											'rot_x', 'rot_x_derivative1', 'rot_x_power2', 'rot_x_derivative1_power2', 
											'rot_y', 'rot_y_derivative1', 'rot_y_power2', 'rot_y_derivative1_power2', 
											'rot_z', 'rot_z_derivative1', 'rot_z_derivative1_power2', 'rot_z_power2']]

				csf = confounds['csf']
				whiteMatter = confounds['white_matter']
				globalSignal = confounds['global_signal']
				dvars = confounds['dvars']
				framewiseDisplacement = confounds['framewise_displacement']
				motionOutlier = confounds['motion_outlier00']

				# build list of confounds to test
				confoundsToTest.append(np.nan_to_num(np.array(sixMotion)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(sixCombinedACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))
				confoundsToTest.append(np.nan_to_num(np.array(dvars)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(sixCombinedACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))
				confoundsToTest.append(np.nan_to_num(np.array(framewiseDisplacement)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(sixCombinedACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))

				confoundsToTestNames.append('sixMotion')
				confoundsToTestNames.append('sixMotion-aCompCor00')
				confoundsToTestNames.append('sixMotion-tCompCor00')
				confoundsToTestNames.append('sixMotion-csf-wm')
				confoundsToTestNames.append('dvars')
				confoundsToTestNames.append('dvars-aCompCor00')
				confoundsToTestNames.append('dvars-tCompCor00')
				confoundsToTestNames.append('dvars-csf-wm')
				confoundsToTestNames.append('framewiseDisplacement')
				confoundsToTestNames.append('framewiseDisplacement-aCompCor00')
				confoundsToTestNames.append('framewiseDisplacement-tCompCor00')
				confoundsToTestNames.append('framewiseDisplacement-csf-wm')

			else:
				sixMotion = confounds[['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']]
				csf = confounds['csf']
				whiteMatter = confounds['white_matter']
				globalSignal = confounds['global_signal']
				dvars = confounds['dvars']
				framewiseDisplacement = confounds['framewise_displacement']
				sixTCompCor = confounds[['t_comp_cor_00','t_comp_cor_01','t_comp_cor_02','t_comp_cor_03','t_comp_cor_04','t_comp_cor_05']]
				sixACompCor = confounds[['a_comp_cor_00','a_comp_cor_01','a_comp_cor_02','a_comp_cor_03','a_comp_cor_04','a_comp_cor_05']]

				# build list of confounds to test
				confoundsToTest.append(np.nan_to_num(np.array(sixMotion)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(sixACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(sixMotion)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))
				confoundsToTest.append(np.nan_to_num(np.array(dvars)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(sixACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(dvars)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))
				confoundsToTest.append(np.nan_to_num(np.array(framewiseDisplacement)))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(sixACompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(sixTCompCor.iloc[:,0])))))
				confoundsToTest.append(np.column_stack((np.nan_to_num(np.array(framewiseDisplacement)), np.nan_to_num(np.array(csf)), np.nan_to_num(np.array(whiteMatter)))))

				confoundsToTestNames.append('sixMotion')
				confoundsToTestNames.append('sixMotion-aCompCor00')
				confoundsToTestNames.append('sixMotion-tCompCor00')
				confoundsToTestNames.append('sixMotion-csf-wm')
				confoundsToTestNames.append('dvars')
				confoundsToTestNames.append('dvars-aCompCor00')
				confoundsToTestNames.append('dvars-tCompCor00')
				confoundsToTestNames.append('dvars-csf-wm')
				confoundsToTestNames.append('framewiseDisplacement')
				confoundsToTestNames.append('framewiseDisplacement-aCompCor00')
				confoundsToTestNames.append('framewiseDisplacement-tCompCor00')
				confoundsToTestNames.append('framewiseDisplacement-csf-wm')

			# run fsl_glm
			for conf in range(len(confoundsToTest)):
				confoundMatrix = np.array(confoundsToTest[conf])
				confoundName = confoundsToTestNames[conf]
				outConfoundFilename = eventsFilename.replace("events.tsv","confound-%s.txt"%confoundName)
				outConfoundFile = os.path.join(args.output_dir,"sub-%s"%subjectLabel,"temp",outConfoundFilename)
				np.savetxt(outConfoundFile, confoundMatrix, delimiter='\t')
				confoundev_files = [('confoundev_files(1)','"%s"'%outConfoundFile)]
				design.update(confoundev_files)

				if not os.path.exists(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,"confound-%s"%confoundName)):
					os.makedirs(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,"confound-%s"%confoundName))

				if not os.path.exists(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,"confound-%s"%confoundName,"%s"%runString)):
					os.makedirs(os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,"confound-%s"%confoundName,"%s"%runString))
					
				featSubDir = os.path.join(args.feat_dir,"task-%s"%taskLabel,"sub-%s"%subjectLabel,"confound-%s"%confoundName,"%s"%runString)

				featOutputDir = [('outputdir','"%s"'%featSubDir)]
				design.update(featOutputDir)

				confoundPlotFilename = eventsFilename.replace("events.tsv","confound-%s.svg"%confoundName)
				confoundPlotFile = os.path.join(args.output_dir,"sub-%s"%subjectLabel,"figures",confoundPlotFilename)
				confound_plot(timeScale,confoundMatrix,confoundPlotFile)

				trial = []
				for key,value in design.items():
					trial.append("set fmri({}) {}".format(key,value))
				trial = [sub.replace('set fmri(feat_files(1))','set feat_files(1)') for sub in trial]
				trial = [sub.replace('set fmri(confoundev_files(1))','set confoundev_files(1)') for sub in trial]

				outDesignFilename = outConfoundFilename.replace(".txt","_design.fsf")
				outDesignFile = os.path.join(featSubDir,outDesignFilename)

				with open(outDesignFile, 'w') as f:
					for item in trial:
						f.write("%s\n" % item)

				cmdFEATModel = ("feat_model %s %s"%(os.path.splitext(outDesignFile)[0], outConfoundFile))
				print(cmdFEATModel)
				run(cmdFEATModel)

				designFile = os.path.join(featSubDir,"design.mat")
				fslGLMFilename = os.path.split(funcFile)[-1].replace("desc-preproc_bold.nii.gz","mask-%d_confound-%s_desc-PEs_stats.nii.gz"%(roiRegionNumber,confoundName))
				fslGLMFile = os.path.join(featSubDir,fslGLMFilename)

				cmdFSLGLM = ("fsl_glm -i %s -d %s --demean -m %s -o %s"%(funcFile,designFile,roiResampleFile,fslGLMFile))
				print(cmdFSLGLM)
				run(cmdFSLGLM)

				glmFile = load_img(fslGLMFile)
				peData = glmFile.get_data()
				peData2D = np.reshape(peData,((np.shape(peData)[0]*np.shape(peData)[1]*np.shape(peData)[2]),np.shape(peData)[3]))
				peData2Dnonzero = []
				for pe in range(np.shape(peData2D)[1]):
					peData2DTemp = peData2D[:,pe]
					peData2Dnonzero.append(peData2DTemp[np.nonzero(peData2DTemp)])
				peData2Dnonzero = np.transpose(np.array(peData2Dnonzero))
				rssTemp = np.zeros(np.shape(peData2Dnonzero))
				tssTemp = np.zeros(np.shape(peData2Dnonzero))

				for pe in range(np.shape(peData2Dnonzero)[1]):
					for voxel in range(np.shape(peData2Dnonzero)[0]):
						rssTemp[voxel,pe] = np.square(peData2Dnonzero[voxel,pe] - np.mean(peData2Dnonzero))
						tssTemp[voxel,pe] = np.square(peData2Dnonzero[voxel,pe])
				rss = np.sum(np.sum(rssTemp,axis=0),axis=0)
				tss = np.sum(np.sum(tssTemp,axis=0),axis=0)
				rssRun.append(rss)
				tssRun.append(tss)

	rssRun = np.array(rssRun)
	rssRun.reshape((-1,len(confoundsToTest)))
	tssRun = np.array(tssRun)
	tssRun.reshape((-1,len(confoundsToTest)))
	outRSSFilename = eventsFilename.replace("events.tsv","RSS.txt"%confoundName)
	outRSSFile = os.path.join(args.output_dir,"sub-%s"%subjectLabel,outRSSFilename)
	outTSSFilename = eventsFilename.replace("events.tsv","TSS.txt"%confoundName)
	outTSSFile = os.path.join(args.output_dir,"sub-%s"%subjectLabel,outTSSFilename)
	np.savetxt(outRSSFile,rssRun,delimiter='\t')
	np.savetxt(outTSSFile,tssRun,delimiter='\t')		





	