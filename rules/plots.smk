rule concatenate_jsons:
	input: 
		cosine_dictionaries = lambda wildcards: expand(join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','cosine.json'),confound_name=confound_names,**wildcards)
	output:
		merged_cosine_dictionary = join(confounder_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json')
	script: '../scripts/merge_dictionaries.py'

rule concatenate_mean_jsons:
	input: 
		cosine_dictionaries = lambda wildcards: expand(join(confounder_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','mean_cosine.json'),confound_name=confound_names,**wildcards)
	output:
		merged_cosine_dictionary = join(confounder_dir,'task-{task}',f'{subj_sess_prefix}_mean_cosine.json')
	script: '../scripts/merge_dictionaries.py'

rule plot_cosine:
	input:
		merged_cosine_dictionary = join(confounder_dir,'task-{task}',f'{subj_sess_prefix}_cosine.json')
	output:
		cosine_plot_svg = join(confounder_dir,'figures','task-{task}',f'{subj_sess_prefix}_cosine.svg')
	script: '../scripts/boxplot.py'

rule plot_correlation:
	input:
		bold_signal_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-noconfounds','run-{run}','bold_signal.txt'),
		confounds_tsv = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{{task}}_run-{{run}}_desc-confounds_regressors.tsv'),
		confounds_dictionary = join(confounder_dir,'confounds_dictionary.json')
	params:
		trial_names = config['task_params']['trial_names']
	output:
		correlation_plot_svg = join(confounder_dir,'figures','task-{task}',f'{subj_sess_prefix}_run-{{run}}_correlation.svg')
	script: '../scripts/heatmap_plot.py' 
