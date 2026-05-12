#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 10G
#SBATCH -e slurm.sj.%j.err
#SBATCH -o slurm.sj.%j.out

###---Author : Maria Kondili
###---Date   : 13/01/2026, update:25/02/2026
###---Subject: Run majiq-v3 in group1 -vs- group2 mode 


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


###
### ANNOTATION
###

# usage: majiq-v3 gff3 [-h] [--version] [--license LICENSE] [--features-default | --features-none | --features-tsv tsv] [--types-genes [T ...]] [--types-transcripts [T ...]] [--types-exons [T ...]]
#                      [--types-silent [T ...]] [--types-hard-skip [T ...]] [--overwrite] [--logger LOGGER] [--logfile-only] [--silent] [--debug]
#                      gff3 splicegraph


##----Download REFERENCE FILE::GFF3_ENSEMBL----## 
#curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz -o gencode.v38.annotation.gff3.gz
#> decompress
#gzip -d gencode.v38.annotation.gff3.gz


ensembl_gff="/shared/bank/homo_sapiens/GRCh38.p14/Ensembl_110/gff3/Homo_sapiens.GRCh38.110.gff3"
#or : gtf="/shared/bank/homo_sapiens/GRCh38/gtf/Homo_sapiens.GRCh38.101.gtf"

mkdir -p ${work_dir}/annotations/


##> RUN ANNOTATION: 

#srun pixi run majiq-v3 gff3 --license ${lic} ${ensembl_gff} ${work_dir}/annotations/sg_hg38.zarr 


###
###  SJ (splice-junctions)
###


##USAGE: srun pixi run majiq-v3 sj [-h] [--version] [--license LICENSE] [--prefix PREFIX] [--no-update-exons] [--update-minreads N] [--update-minpos N] [--strandness {AUTO,NONE,FORWARD,REVERSE}]
                   # [--auto-minreads N] [--auto-minjunctions N] [--auto-mediantolerance X] [--allow-disjoint-contigs | --reject-disjoint-contigs] [--nthreads N]
                   # [--work-directory path | --tmp-work-directory-parent path] [--overwrite] [--logger LOGGER] [--logfile-only] [--silent] [--debug]
                   # bam splicegraph sj

##> OUTPUT: folder sg_hg38.zarr >> contigs, exons,genes,etc..


bam_dir="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/my_STAR_pipeline/Outputs/MarkDuplicates/"
mkdir -p ${work_dir}/SJs/

srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3K_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3K.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3L_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3L.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3M_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3M.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3N_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3N.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3P_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3P.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3Q_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3Q.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3R_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3R.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3S_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3S.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3U_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3U.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3W_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3W.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3Y_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3Y.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I3Z_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I3Z.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I44_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I44.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I45_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I45.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I47_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I47.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I4F_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I4F.sj
srun pixi run majiq-v3 sj --license ${lic} --strandness REVERSE --nthreads ${SLURM_CPUS_PER_TASK} ${bam_dir}/C002I4G_Aligned.sortedByCoord.mrkDups.bam ${work_dir}/annotations/sg_hg38.zarr ${work_dir}/SJs/C002I4G.sj


