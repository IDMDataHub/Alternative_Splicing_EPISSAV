#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --array 0-16
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH -e slurm.build.%j.err
#SBATCH -o slurm.build.%j.out

###---Author : Maria Kondili
###---Date   : 05/03/2026
###---Subject: Run majiq-v3 in group1 -vs- group2 mode for all 16 samples (K ran in test)


###---> Docs :
### https://biociphers.bitbucket.io/majiq-docs/commandbuilder.html
### https://biociphers.bitbucket.io/majiq-docs/getting-started-guide/quantifiers.html
### https://biociphers.bitbucket.io/majiq-docs/v2-to-v3.html

module load pixi
# go in the pixi workspace
work_dir="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_wsp/"
cd ${work_dir}

## Run inside pixi-env : $ pixi shell


###
### LICENSE : 
###
lic="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_license_academic_official.lic"


## Create the arrays for each sample: 
SAMPLES=(C002I3K C002I3L C002I3M C002I3N C002I3P C002I3Q C002I3R C002I3S C002I3U C002I3W C002I3Y C002I3Z C002I44 C002I45 C002I47 C002I4F C002I4G)

CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

###
### SPLICEGRAPH per Sample (created already since 1st run of majiq for each sample)
###


# mkdir -p ${work_dir}/SJs/

# srun pixi run majiq-v3 sj --license ${lic} \
# 							--strandness REVERSE \
# 							--nthreads 4 ${bam_dir}/${SAMPLES[$SLURM_ARRAY_TASK_ID]}_Aligned.sortedByCoord.mrkDups.bam \
# 							${work_dir}/annotations/sg_hg38.zarr \
# 							${work_dir}/SJs/${SAMPLES[$SLURM_ARRAY_TASK_ID]}.sj



###
### BUILD 
###



# the main build command (similar to v2 majiq build)
# the config file is specified as a TSV with a 'group', 'prefix' (experiment name), and 'sj' (path to sj file) columns

# srun pixi run majiq-v3 build 

#> OUTPUT folder .zarr : contigs/, exons/, genes/, introns/, junctions/, .zgroup , .zmetadata


mkdir -p ${work_dir}/build_output/



srun pixi run majiq-v3 build --license ${lic} \
						--nthreads ${SLURM_CPUS_PER_TASK}  \
 						${work_dir}/annotations/sg_hg38.zarr \
 						${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
 						--groups-tsv ${work_dir}/group_configs/config_${CURRENT_SAMPLE}.tsv \
 						 --overwrite
