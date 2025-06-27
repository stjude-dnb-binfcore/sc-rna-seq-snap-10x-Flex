#!/bin/bash

set -e
set -o pipefail

########################################################################
# Load modules
#module load python/3.9.9
module load cellranger/8.0.1

# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Read root path
rootdir=$(realpath "./../..")
echo "$rootdir"

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
cellranger_parameters=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'cellranger_parameters:' | awk '{print $2}')
cellranger_parameters=${cellranger_parameters//\"/}  # Removes all double quotes
echo "$cellranger_parameters"  # Output: This is a string with quotes.

multi_config_dir=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'multi_config_dir:' | awk '{print $2}')
multi_config_dir=${multi_config_dir//\"/}  # Removes all double quotes
echo "multi_config dir: $multi_config_dir"  # Output 

multi_config_file=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'multi_config_file:' | awk '{print $2}')
multi_config_file=${multi_config_file//\"/}  # Removes all double quotes
echo "multi_config file: $multi_config_file"  # Output 


########################################################################
# Create output directories
module_dir="${rootdir}/analyses/cellranger-analysis"
results_dir="${module_dir}/results"
logs_dir="${results_dir}/01_logs"
output_dir="${results_dir}/02_cellranger_count/${cellranger_parameters}"

mkdir -p "$logs_dir" "$output_dir"

echo "Module dir: $module_dir"
echo "Results dir: $results_dir"
echo "Logs dir: $logs_dir"

########################################################################
# multi_config CSV (used by Cell Ranger Multi)
SAMPLES_FILE="${multi_config_dir}/${multi_config_file}"

# Validate sample file
if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Error: CSV file not found at $SAMPLES_FILE"
  exit 1
fi

########################################################################
# STEP 1 ###############################################################
########################################################################
# Submit Cell Ranger Multi job

run_id="multi_run_${cellranger_parameters}"

# Minimum System Requirements based on 
# https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements

bsub -J "RNA_Multi" -n 8 -M 128000 -R "rusage[mem=16000]" \
  -o "${logs_dir}/${run_id}.out" -e "${logs_dir}/${run_id}.err" \
  "cd ${output_dir} && \
   cellranger multi \
     --id=${run_id} \
     --csv=${SAMPLES_FILE} \
     --localcores=8 \
     --localmem=128 \
     --jobmode=lsf"

if [ $? -eq 0 ]; then
  echo "Cell Ranger Multi job submitted successfully."
else
  echo "Error submitting Cell Ranger Multi job."
fi

################################
