run plot_cosine:
	input:
		cosine_dictionary = join(confounder_dir,subj_sess_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json')
	output:
		cosine_plot_svg = join(confounder_dir,subj_sess_dir,figures,'task-{task}',f'{subj_sess_prefix}_cosine.svg')
	script:
		boxplot {input.cosine_dictionary} {output.cosine_plot_svg}

run plot_correlation:
	input:
		bold_signal_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}',f'bold_signal.txt')
		confounds_tsv = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{task}_run-{run}_desc-confounds_regressors.tsv')
	output:
		correlation_plot_svg = join(confounder_dir,subj_sess_dir,figures,'task-{task}',f'{subj_sess_prefix}_run-{run}_correlation.svg')
	params:
		confounds = config['CONFOUNDS']
		trial_names = config['task_params']['trial_names']
	script:
		heatmap_plot {input.confounds_tsv} {input.bold_signal_file} {output.correlation_plot_svg} --confounds_use {params.confounds} --trial_names {params.trial_names}
