if not config["MASK_FILE"]:
	rule run_fsl_glm:
		input:
			func_file = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{{task}}_run-{{run}}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'),
			design_mat = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.mat')
		output:
			stats_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),
			bold_signal_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','bold_signal.txt')
		shell:
			"fsl_glm -i {input.func_file} -d {input.design_mat} --demean -o {output.stats_file} &&"
			"Vest2Text {input.design_mat} {output.bold_signal_file}"
else:
	rule run_fsl_glm_mask:
		input:
			func_file = join(fmriprep_dir,subj_sess_dir,'func',f'{subj_sess_prefix}_task-{{task}}_run-{{run}}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'),
			design_mat = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','design.mat'),
			mask_file = config["MASK_FILE"]
		output:
			stats_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','stats.nii.gz'),
			bold_signal_file = join(feat_dir,'task-{task}',subj_sess_dir,'confound-{confound_name}','run-{run}','bold_signal.txt')
		shell:
			"fsl_glm -i {input.func_file} -d {input.design_mat} --demean -m {input.mask_file} -o {output.stats_file} &&"
			"Vest2Text {input.design_mat} {output.bold_signal_file}"