import os
import fsl.data.featanalysis as FA

# get full path to $FSLDIR
FSLDIR = os.path.expandvars('$FSLDIR')

# get path of output design.fsf file
feat_subdir = os.path.split(snakemake.output.output_design_fsf)[0]

# set variables

# load in design.fsf
design = FA.loadSettings(snakemake.input.feat_dir)

for ev in range(len(snakemake.input.event_files)):
    evTemp = ev + 1
    evName = design['evtitle%d'%evTemp]
    evFile = [evf for evf in snakemake.input.event_files if evf == 'trial-%s_onsets.txt'%evName]
    custom = [('custom%d'%evTemp,'"%s"'%evFile)]
    design.update(custom)

# begin tedious process of rebuilding design.fsf 
feat_files = [('feat_files(1)','"%s"'%snakemake.input.func_file)]
design.update(feat_files)
npts = [('npts','%d'%snakemake.params.bold_reps)]
design.update(npts)
tr = [('tr','%f'%snakemake.params.tr)]
design.update(tr)

feat_output_dir = [('outputdir','"%s"'%feat_subdir)]
design.update(feat_output_dir)

featwatcher = [('featwatcher_yn',0)]
design.update(featwatcher)

regstandard = [('regstandard','"%s/data/standard/MNI152_T1_2mm_brain"'%FSLDIR)]
design.update(regstandard)

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

# if confound_name is 'none' skip step to add confound_file
if snakemake.params.confound_name == 'noconfounds':
    trial = []
    for key,value in design.items():
	    trial.append("set fmri({}) {}".format(key,value))
	    trial = [sub.replace('set fmri(feat_files(1))','set feat_files(1)') for sub in trial]

    with open(snakemake.output.design_fsf, 'w') as f:
	    for item in trial:
	        f.write("%s\n" % item)
else:
    confoundevs = [('confoundevs',1)]
    design.update(confoundevs)
    confoundev_files = [('confoundev_files(1)','"%s"'%snakemake.input.confound_file)]
    design.update(confoundev_files)
    
    trial = []
    for key,value in design.items():
	    trial.append("set fmri({}) {}".format(key,value))
	    trial = [sub.replace('set fmri(feat_files(1))','set feat_files(1)') for sub in trial]

    with open(snakemake.output.design_fsf, 'w') as f:
	    for item in trial:
	        f.write("%s\n" % item)