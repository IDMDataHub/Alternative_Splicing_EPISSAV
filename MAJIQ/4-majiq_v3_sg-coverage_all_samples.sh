#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --array 0-16
#SBATCH --cpus-per-task 6
#SBATCH --mem 10G
#SBATCH -e slurm.sg-cov.%j.err
#SBATCH -o slurm.sg-cov.%j.out

###---Author : Maria Kondili
###---Date   : 24/03/2026
###---Subject: Run sg-coverage for all samples, needed for voila-view & modulizer to detect splicing-categories


module load pixi
# go in the pixi workspace
work_dir="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_wsp"
cd ${work_dir}


lic="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_license_academic_official.lic"


SAMPLES=(C002I3K C002I3L C002I3M C002I3N C002I3P C002I3Q C002I3R C002I3S C002I3U C002I3W C002I3Y C002I3Z C002I44 C002I45 C002I47 C002I4F C002I4G)


CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Current Sample: ${CURRENT_SAMPLE}"

##-- other samples to be -vs-CURRENT--
SJ_FILES=();
sampleIDs=();
for sample in "${SAMPLES[@]}"; do
    if [[ "$sample" != "$CURRENT_SAMPLE" ]]; then
        SJ_FILES+=("${work_dir}/SJs/${sample}.sj")
        sampleIDs+=("${sample}")
    fi
done


echo "Other SampleIDs:"
echo ${sampleIDs[@]}


###
### SG-COVERAGE [.sgc]
###


mkdir -p ${work_dir}/sg_coverage/

# usage: majiq-v3 sg-coverage [-h] [--version] [--license LICENSE]
#                             [--prefixes PREFIXES [PREFIXES ...]]
#                             [--chunksize CHUNKS] [--nthreads N]
#                             [--work-directory path | --tmp-work-directory-parent path]
#                             [--overwrite] [--logger LOGGER] [--logfile-only]
#                             [--silent] [--debug]
#                             splicegraph sg_coverage sj [sj ...]


# for principal-patient
srun pixi run majiq-v3 sg-coverage  --license ${lic} \
									--nthreads ${SLURM_CPUS_PER_TASK} \
                                    ${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
                                    ${work_dir}/sg_coverage/${CURRENT_SAMPLE}_vs_all.sgc \
                                    ${work_dir}/SJs/${CURRENT_SAMPLE}.sj \
                                     

#for all other samples "vs-principal-patient"
srun pixi run majiq-v3 sg-coverage  --license ${lic} \
									--nthreads ${SLURM_CPUS_PER_TASK} \
                                    ${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
                                    ${work_dir}/sg_coverage/all_vs_${CURRENT_SAMPLE}.sgc  \
                                    "${SJ_FILES[@]}" \
                                    


###
### VIEW
###


##--> To run in local Server, or Ask IFB for permission to access host 127.0.0.1:8000
## instructions: https://biociphers.bitbucket.io/majiq/VOILA_view_server.html

# usage: voila view [-h] [--ignore-inconsistent-group-errors] [--long-read-file LONG_READ_FILE] [--clin-controls-file CLIN_CONTROLS_FILE]
#                   [--psicov-grouping-file PSICOV_GROUPING_FILE] [--group-order-override-file GROUP_ORDER_OVERRIDE] 
#					[-p PORT] [--host HOST]
#                   [--web-server {waitress,gunicorn,flask}] [--num-web-workers NUM_WEB_WORKERS] [--gunicorn-worker-class GUNICORN_WORKER_CLASS] [--enable-passcode]
#                   [--index-file INDEX_FILE] [--force-index] [--enable-type-indexing] [--strict-indexing] [--only-index] [-j NPROC] [--debug] [--lazy-load-zarr]
#                   [--parallel-chunksize PARALLEL_CHUNKSIZE] [-l LOGGER] [--silent]
#                   files [files ...]


# srun pixi run voila view -p 8000 --host 127.0.0.1 --web-server gunicorn  -j ${SLURM_CPUS_PER_TASK}  \
#                           ${work_dir}/build_output/sg_${CURRENT_S[$SLURM_ARRAY_TASK_ID]}_vs_all.zarr \
#                           ${work_dir}/sg_coverage/${CURRENT_SAMPLE}.sgc \
#                           ${work_dir}/sg_coverage/all-vs-${CURRENT_S[$SLURM_ARRAY_TASK_ID]}.sgc \
#                           ${work_dir}/heterogen/patientK_vs_controls.heterogen.voila
