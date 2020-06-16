from os.path import join
import json
from snakemake.utils import format
from bids import BIDSLayout
import pandas as pd
import numpy as np

configfile: "config.yaml"

if not os.path.isdir(config["BIDS_DIR"]):
    raise Exception("Cannot find the dataset at %s" % config["BIDS_DIR"])

# can also override at command-line with e.g.:  --config bids_dir='path/to/dir'  or --configfile ...
bids_dir = config['BIDS_DIR']
feat_dir = config['FEAT_DIR']
fmriprep_dir = join(config['FMRIPREP_DIR'],'fmriprep')
confounder_dir = config['CONFOUNDER_DIR']

# pull confounds dictionary from config file and save
# create list of confound grouping names and prepend
# with "noconfounds" to account for testing GLM with
# no experimental confounds
confounds = config['CONFOUNDS']
confound_names = list(confounds.keys())
confound_names.insert(0,'noconfounds')

# save out config file CONFOUNDS dictionary as json file
confounds_dictionary = join(confounder_dir,'confounds_dictionary.json')
with open(confounds_dictionary, 'w') as fp:
	json.dump(confounds, fp, indent=4, separators=(',', ': '))

# load participants.tsv file, and strip of sub- from participant_id column
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.str.slice(4).to_list() 

layout = BIDSLayout(fmriprep_dir,validate=False)

# get entities from cfg file to select input files from fmriprep
fmriprep_params = config['fmriprep_params']
sessions = layout.get_sessions()
#runs_list = layout.get_runs()
runs_list = config['fmriprep_params']['num_runs']
runs_list = np.arange(runs_list)+1

# convert num_runs into list of strings matching "run-%02d" in BIDS filenames
runs_list = np.array(runs_list)
runs = []
for r in range(len(runs_list)):
    x = runs_list[r]
    y = (np.array2string(x, formatter={'int':lambda x: "%02d" % x}))
    runs.append(y)

# get task parameters
task = config['fmriprep_params']['task']
trial_names = config['task_params']['trial_names']

# create strings including wildcards for subj_sess_dir and subj_sess_prefix
if len(sessions) > 0:
    subj_sess_dir = join('sub-{subject}','ses-{session}')
    subj_sess_prefix = 'sub-{subject}_ses-{session}'
else:
    subj_sess_dir = 'sub-{subject}'
    subj_sess_prefix = 'sub-{subject}'

rule all:
    input: expand(join(confounder_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json'),task=task,subject=subjects,session=sessions,allow_missing=True)
        
include: 'rules/fsl_design.smk'
include: 'rules/fsl_glm.smk'
include: 'rules/cosine_angle.smk'
include: 'rules/plots.smk'
