#!/usr/bin/env bash
#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH -p fast
#SBATCH -o slurm.rmats.%j.out
#SBATCH -e slurm.rmats.%j.err

#module load conda
module load python/3.9
module load rmats/4.3.0 #> since 2024
#module load r/4.4.1 

rmats_dir="/shared/ifbstor1/software/miniconda/envs/rmats-4.3.0/bin"
cores=${SLURM_CPUS_PER_TASK}
gtf="/shared/bank/homo_sapiens/GRCh38.p14/Ensembl_110/gtf/Homo_sapiens.GRCh38.110.gtf" #> same as the STAR & RSEM 
workDir="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/rMATS/"


##--> Input files, given in cmd-arguments: 
##--> for each pair of sample/all.wo.sample, run same script in a loop)
all_samples="${workDir}/input_bam/all_samples.csv"


##--> Create output-directory 
outDir="${workDir}/Outputs/one_group_stat_off/output_all_statoff"; 
mkdir -p ${outDir}
tmpDir="${outDir}/tmp_all_statoff/"; 
mkdir -p ${tmpDir}


##--> Run rmats --STATOFF for all in ONE GROUP

## python ${rmats_dir}/rmats.py \
rmats.py --b1 ${all_samples} \
--gtf ${gtf} --readLength 101 \
--od ${outDir} --tmp ${tmpDir} \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} \
--statoff \
--task both ## ? why both if stat=off ? 

##>30/09,11:22 , JOBID = 61727842

##--> Run STAT model separetely 

# rmats.py --od ${outDir} \
# 		--tmp ${tmpDir} \
# 		--task stat 


### ATTENTION TO Parameters :
### --libType = fr-firststrand = for "reverse" strand, 
###           = fr-secondstrand = for "forward" strand
##>> About stranded-RNA protocol : https://chipster.csc.fi/manual/library-type-summary.html 
### --readLength = 100 -> seen in FastQC-report (x-axis of quality-per-bp)


## run with :
##--> sbatch run_rmats_both.sh $1=input_bam/sample_K.csv $2=input_bam/samples_K_vs_all.csv

