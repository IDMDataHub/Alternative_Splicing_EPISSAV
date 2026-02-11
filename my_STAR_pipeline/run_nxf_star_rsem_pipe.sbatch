#!/usr/bin/env bash
#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --cpus-per-task=20
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load nextflow/24.04.4

work_dir=/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/my_STAR_pipeline/
script_dir=${work_dir}/NXF_workspace  #< where main.nf is
nxf_wsp=${script_dir}/Pipe_Output

# Create if it doesn't exist:
mkdir -p ${nxf_wsp}

# cd ${work_dir}/Nextflow_Pip/

###>> Change of Command to run Nextflow since 2021
#1 / config file must be in same dir as .nf file ,named "nextflow.config"
# 2/ .nf script should be named "main.nf"
###  mkdir rnaseq_pip
###  mv config_pip.conf rnaseq_WSP/nextflow.config
###  mv rnaseq_nxf_pip.nf  rnaseq_WSP/main.nf
# 3/ above files must be in a dir named after the project or Pipeline ,for this one : rnaseq_WSP/

##increase the max number of open files with Linux command ulimit -n
# ulimit -n 200000

##> call Nextflow with -w WORKSPACE SCRIPTS_DIR with relative path,not root-path
##> 17/09 : with MarkDuplicates-process added
nextflow  run -w ${nxf_wsp} ${script_dir} --samplelist ${work_dir}/Samplesheets/all_patients.tsv -resume

