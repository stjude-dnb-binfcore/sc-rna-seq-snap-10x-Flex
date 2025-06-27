# Pipeline for summarizing results from all modules and generating project reports for sc-/sn-RNA-Seq Analysis in 10X Flex Genomics data

## Usage

`run-project-updates.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.

The user should edit the last section of the report, labeled `# Future Directions`, in the `01-generate-project-report.Rmd` script after reviewing the results and adding their recommendations for next steps. 

Links are provided in each analysis section. If the module has been run and results are available in the folder, the link will work. However, if the analysis module hasn't been run, attempting to access the link will result in an error.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node (while within the container):

```
bash run-project-updates.sh
```

### Run module by using lsf on HPC with the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content

This folder contains scripts tasked for summarizing results from all modules and generating project reports for sc-/sn-RNA-Seq Analysis in 10X Genomics data.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-generate-project-report.Rmd
├── lsf-script.txt
├── README.md
├── reports
├── run-project-updates.R
└── run-project-updates.sh
```
