if not config['CONTRAST']:
	rule run_cosine_angle:
		input: 
			stats_files = lambda wildcards: expand(join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),run=runs,**wildcards)
		params:
			trial_names = config['task_params']['trial_names']
		output:
			cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}',f'{subj_sess_prefix}_cosine.json'),
			mean_cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}',f'{subj_sess_prefix}_mean_cosine.json')
		script: '../scripts/cosine_angle.py'       	
else:
	rule run_cosine_angle_contrast:
		input: 
			stats_files = lambda wildcards: expand(join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),run=runs,**wildcards)
		params:
			trial_names = config['task_params']['trial_names'],
			contrast_vector = config['CONTRAST']
		output:
			cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}',f'{subj_sess_prefix}_cosine.json'),
			mean_cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}',f'{subj_sess_prefix}_mean_cosine.json')
		script: '../scripts/cosine_angle_contrast.py'
	