#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --array 0-16
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH -e slurm.psicov.%j.err
#SBATCH -o slurm.psicov.%j.out

###---Author : Maria Kondili
###---Date   : 19/03/2026
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


##--- Create the array for the samples: 
SAMPLES=(C002I3K C002I3L C002I3M C002I3N C002I3P C002I3Q C002I3R C002I3S C002I3U C002I3W C002I3Y C002I3Z C002I44 C002I45 C002I47 C002I4F C002I4G)


##----Create the array of the SJs excluding the current sample, to compare them to it on the command 

CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Current Sample: ${CURRENT_SAMPLE}"

SJ_FILES=()
sampleIDs=()
for sample in "${SAMPLES[@]}"; do
    if [[ "$sample" != "$CURRENT_SAMPLE" ]]; then
        SJ_FILES+=("SJs/${sample}.sj")
        sampleIDs+=("${sample}")
    fi
done

echo "Other SampleIDs:"
echo ${sampleIDs[@]}

###
###  PSI COVERAGE
###
##--> general command by array


# usage: majiq-v3 psi-coverage [-h] [--version] [--license LICENSE]
#                              [--prefixes PREFIXES [PREFIXES ...]]
#                              [--minreads QUANTIFY_MINREADS]
#                              [--minbins QUANTIFY_MINBINS]
#                              [--stack-pvalue-threshold P] [--ignore-from sg]
#                              [--strict-lsvs | --permissive-lsvs | --source-lsvs | --target-lsvs]
#                              [--nthreads N]
#                              [--work-directory path | --tmp-work-directory-parent path]
#                              [--overwrite] [--logger LOGGER] [--logfile-only]
#                              [--silent] [--debug]
#                              splicegraph psi_coverage sj [sj ...]


mkdir -p ${work_dir}/psi_cov/


if [ -e ${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr ] ; then 


	##example group-1: majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
 	srun pixi run majiq-v3 psi-coverage --license ${lic} --nthreads ${SLURM_CPUS_PER_TASK} \
 										${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
 										${work_dir}/psi_cov/${CURRENT_SAMPLE}_vs_all.psicov \
 										${work_dir}/SJs/${CURRENT_SAMPLE}.sj --overwrite


 	#example group-2: majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Brain_Cerebellum.psicov /path/to/results/sj/sample_bc_1.sj /path/to/results/sj/sample_bc_2.sj
 	srun pixi run majiq-v3 psi-coverage --license ${lic} --nthreads ${SLURM_CPUS_PER_TASK} \
 											${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
 											${work_dir}/psi_cov/all_vs_${CURRENT_SAMPLE}.psicov  \
 											"${SJ_FILES[@]}"  --overwrite

fi 


#---------------------------------------simpler method with manually def. each sample -----------------------------------------------------------------------------------------------------#
##--> M
# if [ -e ${work_dir}/build_output/sg_C002I3M_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3M_vs_all.zarr ${work_dir}/psi_cov/all_vs_M.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 


#  ##--> N

#  if [ -e ${work_dir}/build_output/sg_C002I3N_vs_all.zarr/ ] ; then 
 
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3N_vs_all.zarr ${work_dir}/psi_cov/all_vs_N.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 


#  ##--> P

#  if [ -e ${work_dir}/build_output/sg_C002I3P_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3P_vs_all.zarr ${work_dir}/psi_cov/all_vs_P.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 


#  ##--> Q

#  if [ -e ${work_dir}/build_output/sg_C002I3Q_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3Q_vs_all.zarr ${work_dir}/psi_cov/all_vs_Q.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 


#  ##--> R

#  if [ -e ${work_dir}/build_output/sg_C002I3R_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3R_vs_all.zarr ${work_dir}/psi_cov/all_vs_R.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> S
 
#  if [ -e ${work_dir}/build_output/sg_C002I3S_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3S_vs_all.zarr ${work_dir}/psi_cov/all_vs_S.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> U
 
#  if [ -e ${work_dir}/build_output/sg_C002I3U_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3U_vs_all.zarr ${work_dir}/psi_cov/all_vs_U.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> W
 
#  if [ -e ${work_dir}/build_output/sg_C002I3W_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3W_vs_all.zarr ${work_dir}/psi_cov/all_vs_W.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> Y
 
#  if [ -e ${work_dir}/build_output/sg_C002I3Y_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3Y_vs_all.zarr ${work_dir}/psi_cov/all_vs_Y.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 


#  ##--> Z
 
#  if [ -e ${work_dir}/build_output/sg_C002I3Z_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I3Z_vs_all.zarr ${work_dir}/psi_cov/all_vs_Z.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> 44
 
#  if [ -e ${work_dir}/build_output/sg_C002I44_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I44_vs_all.zarr ${work_dir}/psi_cov/all_vs_44.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> 45
 
# if [ -e ${work_dir}/build_output/sg_C002I45_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I45_vs_all.zarr ${work_dir}/psi_cov/all_vs_45.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

 
#  ##--> 47
 
#  if [ -e ${work_dir}/build_output/sg_C002I47_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I47_vs_all.zarr ${work_dir}/psi_cov/all_vs_47.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I4F.sj \
#  		${work_dir}/SJs/C002I4G.sj --overwrite

#  fi 

#  ##--> 4G

#  if [ -e ${work_dir}/build_output/sg_C002I4G_vs_all.zarr/ ] ; then 
 	
# 	##example : majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj
#  	srun pixi run majiq-v3 psi-coverage --license ${lic} ${work_dir}/build_output/sg_C002I4G_vs_all.zarr ${work_dir}/psi_cov/all_vs_4G.psicov ${work_dir}/SJs/C002I3K.sj \
#  		${work_dir}/SJs/C002I3L.sj \
#  		${work_dir}/SJs/C002I3M.sj \
#  		${work_dir}/SJs/C002I3N.sj \
#  		${work_dir}/SJs/C002I3P.sj \
#  		${work_dir}/SJs/C002I3Q.sj \
#  		${work_dir}/SJs/C002I3R.sj \
#  		${work_dir}/SJs/C002I3S.sj \
#  		${work_dir}/SJs/C002I3U.sj \
#  		${work_dir}/SJs/C002I3W.sj \
#  		${work_dir}/SJs/C002I3Y.sj \
#  		${work_dir}/SJs/C002I3Z.sj \
#  		${work_dir}/SJs/C002I44.sj \
#  		${work_dir}/SJs/C002I45.sj \
#  		${work_dir}/SJs/C002I47.sj \
#  		${work_dir}/SJs/C002I4F.sj --overwrite

#  fi 