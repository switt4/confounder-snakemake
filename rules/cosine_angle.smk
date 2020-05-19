if config['CONSTRAST']:
	rule run_cosine_angle:
		input: 
			stats_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz')
        	param_estimate_files = join(feat_dir,'param_estimate_file_list.txt')
		params:
			trial_names = config['task_params']['trial_names']
			contrast_vector = config['CONTRAST']
			confound_names = confound_names
		output:
			cosine_dictionary = join(confounder_dir,subj_sess_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json')
			mean_cosine_dictionary = join(confounder_dir,subj_sess_dir,'task-{task}',f'{subj_sess_prefix}_mean_cosine.json')
		script: 'scripts/cosine_angle_contrast.py'
        	
else:
	rule run_cosine_angle:
		input: 
			stats_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz')
        	param_estimate_files = join(feat_dir,'param_estimate_file_list.txt')
		params:
			trial_names = config['task_params']['trial_names']
			confound_names = confound_names
		output:
			cosine_dictionary = join(confounder_dir,subj_sess_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json')
			mean_cosine_dictionary = join(confounder_dir,subj_sess_dir,'task-{task}',f'{subj_sess_prefix}_mean_cosine.json')
		script: 'scripts/cosine_angle.py'