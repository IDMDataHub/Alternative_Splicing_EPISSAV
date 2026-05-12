#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --array 0-16
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH -e slurm.deltapsi.%j.err
#SBATCH -o slurm.deltapsi.%j.out

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


SAMPLES=(C002I3K C002I3L C002I3M C002I3N C002I3P C002I3Q C002I3R C002I3S C002I3U C002I3W C002I3Y C002I3Z C002I44 C002I45 C002I47 C002I4F C002I4G)

###
### MAJIQ DeltaPSI 
###


# usage: majiq-v3 deltapsi [-h] [--version] [--license LICENSE] -psi1 PSI1
#                          [PSI1 ...] [-psi2 PSI2 [PSI2 ...]] [-n NAME1 NAME2]
#                          [--min-experiments X] [--splicegraph SG]
#                          [--annotated SG] [--output-tsv TSV]
#                          [--output-voila ZARR]
#                          [--select-grp1-prefixes PREFIX [PREFIX ...] |
#                          --drop-grp1-prefixes PREFIX [PREFIX ...]]
#                          [--select-grp2-prefixes PREFIX [PREFIX ...] |
#                          --drop-grp2-prefixes PREFIX [PREFIX ...]]
#                          [--downsample2] [--rng-seed S]
#                          [--allow-prefix-overlap] [--psibins B]
#                          [--changing-threshold CT]
#                          [--nonchanging-threshold NCT]
#                          [--empirical-prior | --default-prior]
#                          [--prior-minreads R] [--prior-minevents N]
#                          [--prior-iter N] [--nthreads N]
#                          [--work-directory path | --tmp-work-directory-parent path]
#                          [--show-progress | --disable-progress]
#                          [--scheduler-address url | --scheduler-file json]
#                          [--overwrite] [--logger LOGGER] [--logfile-only]
#                          [--silent] [--debug]


mkdir -p ${work_dir}/deltapsi/
mkdir -p ${work_dir}/deltapsi/voila/
mkdir -p ${work_dir}/deltapsi/tsv/

##example: majiq-v3 deltapsi --splicegraph /path/to/results/build/sg.zarr --output-voila /path/to/results/dpsi/Brain_Cerebellum-Muscle_Skeletal.dpsicov --output-tsv /path/to/results/dpsi/Brain_Cerebellum-Muscle_Skeletal.tsv -psi1 /path/to/results/psi/Brain_Cerebellum.psicov -psi2 /path/to/results/psi/Muscle_Skeletal.psicov



CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Current Sample: ${CURRENT_SAMPLE}"

srun pixi run majiq-v3 deltapsi --license ${lic} --nthreads ${SLURM_CPUS_PER_TASK} \
								--splicegraph  ${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
				  				--output-voila ${work_dir}/deltapsi/voila/voila_${CURRENT_SAMPLE}.dpsicov \
				  				--output-tsv ${work_dir}/deltapsi/tsv/dpsi_${CURRENT_SAMPLE}.tsv \
				  				-psi1 ${work_dir}/psi_cov/${CURRENT_SAMPLE}_vs_all.psicov \
				  				-psi2 ${work_dir}/psi_cov/all_vs_${CURRENT_SAMPLE}.psicov 




# ##--> M
# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3M_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_M_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_M_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/M_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_M.psicov 


# ##--> N

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3N_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_N_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_N_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/N_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_N.psicov 


# ##--> P

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3P_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_P_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_P_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/P_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_P.psicov 


# ##--> Q

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3Q_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_Q_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_Q_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/Q_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_Q.psicov 

# ##--> R

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3R_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_R_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_R_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/R_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_R.psicov 

# ##--> S

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3S_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_S_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_S_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/S_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_S.psicov


# ##--> U

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3U_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_U_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_U_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/U_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_U.psicov

# ##--> W

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3W_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_W_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_W_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/W_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_W.psicov

# ##--> Y

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3Y_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_Y_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_Y_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/Y_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_Y.psicov

# ##--> Z

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I3Z_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_Z_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_Z_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/Z_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_Z.psicov

# ##--> 44

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I44_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_44_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_44_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/44_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_44.psicov

# ##--> 45

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I45_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_45_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_45_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/45_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_45.psicov

# ##--> 47

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I47_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_47_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_47_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/47_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_47.psicov

# ##--> 4F

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I4F_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_4F_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_4F_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/4F_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_4F.psicov


# ##--> 4G

# srun pixi run majiq-v3 deltapsi --license ${lic} \
# 								--splicegraph  ${work_dir}/build_output/1sample_vs_all/sg_C002I4G_vs_all.zarr \
# 				  				--output-voila ${work_dir}/deltapsi/voila_4G_vs_all.dpsicov \
# 				  				--output-tsv ${work_dir}/deltapsi/dpsi_4G_vs_all.tsv \
# 				  				-psi1 ${work_dir}/psi_cov/4G_vs_all.psicov \
# 				  				-psi2 ${work_dir}/psi_cov/all_vs_4G.psicov