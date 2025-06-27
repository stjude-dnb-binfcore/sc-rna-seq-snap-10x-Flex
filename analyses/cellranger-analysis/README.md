# Pipeline for Cell Multiplexing with cellranger multi for libraries for sc-/sn-RNA-Seq Analysis in 10X Flex Genomics data

## Usage

`submit-multiple-jobs.sh` is designed to be run as if it were called from this module directory, even when called from outside of this directory.

The `submit-multiple-jobs.sh` script is designed to run the following two steps: 
   - Step 1: To run the `j1.sh` script to align single or multiple libraries in parallel, i.e., `run-cellranger-analysis`. 
   - Step 2: To run `j2.sh` to summarize alignment results, i.e., `summarize-cellranger-analysis`. The latter script will be on hold and executed once all libraries are aligned and `j1.sh` is complete. This is been taken care of by the `waiter.sh` script.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `../../project_parameters.Config.yaml`: define the `multi_config_dir`, `multi_config_file`, `cellranger_parameters`. For a list of genome references maintained and supported by the Bioinformatics Core at DNB, please review the [Genome References wiki page](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/wiki/2.-Genome-References). Please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues) to request the path to the reference genome of preference. Otherwise, specify the path to the reference genome of your preference. 


- `j1.sh`: 
  - `--force-cells`: User can add flags as necessary for their analyses and compare alignments with CellRanger, e.g., by using `--force-cells=8000` to constrain the number of cells to be used for alignment. We recommend running by default and, after careful assessment, editing parameters. We have found that the default parameters set up here work well for most cases.

- `j2.sh`: 
   - Format of metrics and statistics from the CellRanger multi algorithm is different than the one from the CellRanger count, so we will need to convert the files. This is being taken care of in the `j2-prepare-CSV.R` script before running the summary. 

- `multi_config_file`:
    - User needs to define the probe set and information for the library to be demultiplexed. We recommend setting `create_bam_value` to `false` to reduce memory usage per project. Moreover, please note that usually there is not much benefit of creating BAM files from the Flex analysis. Similarly, there is really no general utility of aligning the Flex reads to the genome reference, and that's why for Flex data analysis, the [input of transcriptome reference is optional](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#gene-expression).
    
User also need to define `sample_prefix` with the Sample ID used for the samples of the project. Sample IDs should follow a format like: PREFIX001 (e.g., DYE001, ABC-002, XYZ_003). You can specify multiple prefixes if your project uses more than one.

## Run module on HPC

To run all of the scripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bsub < submit-multiple-jobs.sh
```

Please note that this will run the analysis module outside of the container while submitting lsf job on HPC. This is currently the only option of running the `cellranger-analysis` module. By default, we are using `python/3.9.9` and `cellranger/8.0.1` as available on St Jude HPC.


## Folder content
This folder contains scripts tasked to run, demultiplex, and summarize Cell Ranger multi for libraries for sc-/sn-RNA-Seq Analysis in 10X Flex Genomics data across the project. For more information and updates, please see [Cell Ranger support page](https://www.10xgenomics.com/products/flex-gene-expression).

For more information about best approaches and kits for demultiplexing, see [Bioinformatics Tools for Sample Demultiplexing](https://www.10xgenomics.com/analysis-guides/bioinformatics-tools-for-sample-demultiplexing#bioinformatics-tools-for-sample-demultiplexing).

This module uses CellRanger v8.0.1 for the alignment.


## Folder structure 

The structure of this folder is as follows:

```
├── j1.sh
├── j2.sh
├── j2-prepare-CSV.R
├── README.md
├── results
|   ├── 01_logs
|   ├── 02_cellranger_count
|   |   └── {cellranger_parameters}
|   |        └── multi_run_{cellranger_parameters}
|   └── 03_cellranger_count_summary
|   |   └── {cellranger_parameters}
|   |       └── multi_run_{cellranger_parameters}
├── submit-multiple-jobs.sh
├── util
|   └── summarize_cellranger_results.py
└── waiter.sh
```
