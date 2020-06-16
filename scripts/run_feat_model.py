import os
import subprocess

def run(command, env={}):
	merged_env = os.environ
	merged_env.update(env)
	process = subprocess.Popen(command, stdout=subprocess.PIPE,
							   stderr=subprocess.STDOUT, shell=True,
							   env=merged_env)
	while True:
		line = process.stdout.readline()
		line = str(line, 'utf-8')[:-1]
		print(line)
		if line == '' and process.poll() != None:
			break
	if process.returncode != 0:
		raise Exception("Non zero return code: %d"%process.returncode)

design_dir = os.path.dirname(snakemake.input.design_fsf)

if snakemake.params.confound_name == 'noconfounds':
    command = 'feat_model %s/design'%design_dir
    run(command)
else:
    command = 'feat_model %s/design %s'%(design_dir,snakemake.input.confound_file)
    run(command)