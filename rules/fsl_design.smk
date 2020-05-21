rule create_custom_events:
	input:
		events_tsv = join(bids_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{{task}}_run-{{run}}_events.tsv')
	params:
		trial_name = '{trial_name}'
	output:
		event_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','trial-{trial_name}_onsets.txt')
	script: 'scripts/custom_events_file.py'
	
rule create_confound_file:
	input:
		confounds_tsv = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{{task}}_run-{{run}}_desc-confounds_regressors.tsv'),
		confounds_dictionary = join(confounder_dir,'confounds_dictionary.json')
	params:
		confound_name = '{confound_name}'
	log:
		logfile = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','confounds-{confound_name}.log')
	output:
		confound_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','confounds-{confound_name}.txt')
	script: 'scripts/create_confound_files.py'

rule create_fsl_design:
	input:
		feat_dir = config['FEAT_DIR'],
		func_file = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}-task-{{task}}_run-{{run}}_bold.nii.gz'),
		event_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','trial-{trial_name}_onsets.txt'),
		confound_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','confounds-{confound_name}.txt')
	params:
		confound_name = '{confound_name}',
		bold_reps = config['sequence_params']['bold_reps'],
		TR = config['sequence_params']['TR']
	output:
		design_fsf = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.fsf')
	script: 'scripts/create_fsl_design_confounds.py'

rule estimate_fsl_design:
	input:
		design_fsf = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.fsf')
	output:
		design_fsf = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.fsf'),
		design_cov_png = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design_cov.png'),
		design_cov_ppm = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design_cov.ppm'),
		design_con = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.con'),
		design_frf = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.frf'),
		design_mat = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.mat'),
		design_min = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.min'),
		design_png = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.png'),
		design_ppm = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.ppm'),
		design_trg = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.trg')
	shell:'feat_model {input}'
