import numpy as np
import pandas as pd

def read_tsv(inTSV):
	dataTSV = pd.read_table(inTSV)
	dataTSV.head()
	return dataTSV

# load in events.tsv file
event_file = snakemake.input.events_tsv
task_timings = read_tsv(event_file)

# define trial name variable
ev_name = snakemake.params.trial_name

# extract onsets and durations from events.tsv file for each trial type
# save out as text file in 3-column custom format for design.fsf	
timings = task_timings[task_timings.trial_type == ev_name]
timings = timings[['onset','duration']]
timings = np.column_stack((timings,np.ones(len(timings))))
np.savetxt(snakemake.output.event_file,timings,delimiter='\t')