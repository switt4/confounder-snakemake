import os
import fsl.data.featanalysis as FA
import fsl.utils.run as RUN

# get number of time points using fslinfo (in case of non-uniform runs)
info = RUN.runfsl('fslinfo',snakemake.input.func_file)
info_temp = info.split()
info_dct = {info_temp[i]: info_temp[i+1] for i in range(0,len(info_temp),2)}
num_pts = int(info_dct['dim4'])

# get full path to $FSLDIR
FSLDIR = os.path.expandvars('$FSLDIR')

# get path of output design.fsf file
feat_subdir = os.path.split(snakemake.output.design_fsf)[0]

# set variables

# load in design.fsf
design = FA.loadSettings(snakemake.input.feat_dir)

for ev in range(len(snakemake.input.event_files)):
    evTemp = ev + 1
    evName = design['evtitle%d'%evTemp]
    evFile = snakemake.input.event_files[ev]
    #evFile = [evf for evf in snakemake.input.event_files if evf == 'trial-%s_onsets.txt'%evName]
    custom = [('custom%d'%evTemp,'"%s"'%evFile)]
    design.update(custom)

# begin tedious process of rebuilding design.fsf 
feat_files = [('feat_files(1)','"%s"'%snakemake.input.func_file)]
design.update(feat_files)
npts = [('npts','%d'%num_pts)]
design.update(npts)
tr = [('tr','%f'%snakemake.params.TR)]
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
    trial = [sub.replace('set fmri(confoundev_files(1))','set confoundev_files(1)') for sub in trial]

    with open(snakemake.output.design_fsf, 'w') as f:
	    for item in trial:
	        f.write("%s\n" % item)