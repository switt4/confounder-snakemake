if not config['CONTRAST']:
	rule run_cosine_angle:
		input: 
			stats_files = lambda wildcards: expand(join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),run=runs,**wildcards)
		params:
			trial_names = config['task_params']['trial_names'],
			confound_name = '{confound_name}'
		conda:
			'../envs/nilearn.yaml'
		output:
			cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','cosine.json'),
			mean_cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','mean_cosine.json')
		script: '../scripts/cosine_angle.py'
else:
	rule run_cosine_angle_contrast:
		input: 
			stats_files = lambda wildcards: expand(join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),run=runs,**wildcards)
		params:
			trial_names = config['task_params']['trial_names'],
			confound_name = '{confound_name}',
			contrast_vector = config['CONTRAST']
		conda:
			'../envs/nilearn.yaml'
		output:
			cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','cosine.json'),
			mean_cosine_dictionary = join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','mean_cosine.json')
		script: '../scripts/cosine_angle_contrast.py'
	