# confounder-snakemake
confounder is a BIDS app designed to help researchers determine which set of experimental confounds is best suited for improving the sensitivity of standard task-based GLM analyses.

The app computes the cosine similarity of activation patterns across multiple runs of a given task for user-specified sets of experimental confounds.  Currently the app can only support those experimental confounds estimated by a standard fMRIPrep analysis.  Users have the option to supply a co-registered ROI and/or a contrast vector to further refine the results.  Default options estimate cosine similarity of whole-brain activation patterns for a standard effects-of-interest contrast across all task regressors.

Inputs:
- Valid BIDS directory with multiple runs of a task-based fMRI experiment.  At minimum, the `events.tsv` file for each run must have completed `onsets`, `duration` and `trial_type` columns.
- Valid fMRIPrep directory.  App is currently set to use the `desc-preproc.nii.gz` files.
- Directory containing an example FSL `design.fsf` file for the task of interest.  The timings of the design need to be specified using the "3 column custom" format.

Prerequisites:
- conda
- FSL (WIP: suppport for FSL in a container)

Notes:
- App only works for one task at a time.
- All subjects must have the same number of runs of the task.  (Length of individual runs can be variable.)

## Authors
* Suzanne T. Witt @switt4

## Usage

If you use this workflow in a paper, please do not forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the `config.yaml`. Adjust `config.yml` to configure the workflow execution, and `participants.tsv` to specify your subjects.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -np

Execute the workflow locally via

    snakemake --use-conda -p
    
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.
    
### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html
    
Users may find that they need to install two conda packages (`pygraphviz` and `imagemagick`) in order to have the html reports formateed correctly.

    conda install pygraphviz
    conda install -c conda-forge imagemagick

This report can be forwarded to your collaborators.  N.B. All images are embedded in the report, so this file can become large when running large numbers of subjects/runs.

An example can be seen [here](https://github.com/switt4/confounder-snakemake/blob/master/report.html).

To exit the snakemake environment when finished type: `conda deactivate`.



