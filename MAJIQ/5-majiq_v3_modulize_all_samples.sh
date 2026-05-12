#!/usr/bin/env bash

#SBATCH -A splicing_events_detection_rare_diseases
#SBATCH -p fast
#SBATCH --array 0-16
#SBATCH --cpus-per-task 6
#SBATCH --mem 10G
#SBATCH -e slurm.modulize.%j.err
#SBATCH -o slurm.modulize.%j.out

###---Author : Maria Kondili
###---Date   : 14/04/2026
###---Subject: Run "modulize" for all samples (after sg-coverage)


module load pixi
# go in the pixi workspace :should always run MAJIQ from here where env.is !
work_dir="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_wsp"
cd ${work_dir}


lic="/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/MAJIQ/majiq_license_academic_official.lic"


SAMPLES=(C002I3K C002I3L C002I3M C002I3N C002I3P C002I3Q C002I3R C002I3S C002I3U C002I3W C002I3Y C002I3Z C002I44 C002I45 C002I47 C002I4F C002I4G)



CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Current Sample: ${CURRENT_SAMPLE}"


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
### ΜΟDULIZER :https://biociphers.bitbucket.io/majiq-docs/modulizer/quick-start.html
###

#! command from "Quick-Start"  needs to update !

## GOOGLE GROUP : https://groups.google.com/g/majiq_voila?pli=1

# usage: voila modulize [-h] [--overwrite] [--ignore-inconsistent-group-errors]
#                       [--only-binary] [--untrimmed-exons] [--show-all]
#                       [--heatmap-selection {shortest_junction,max_abs_dpsi}]
#                       [--disable-metadata] [--show-read-counts]
#                       [--show-per-sample-psi]
#                       [--psicov-grouping-file PSICOV_GROUPING_FILE]
#                       [--gene-ids [GENE_IDS ...]]
#                       [--debug-num-genes DEBUG_NUM_GENES] [--output-mpe]
#                       [--putative-multi-gene-regions]
#                       [--keep-constitutive [KEEP_CONSTITUTIVE]] --> default = 1 read to keep constitutive
#                       [--keep-no-lsvs-modules] [--keep-no-lsvs-junctions]
#                       [--decomplexify-psi-threshold DECOMPLEXIFY_PSI_THRESHOLD]			--> default = "0.05"
#                       [--decomplexify-deltapsi-threshold DECOMPLEXIFY_DELTAPSI_THRESHOLD] --> default = "0.0"
#                       [--decomplexify-reads-threshold DECOMPLEXIFY_READS_THRESHOLD] 		--> The default = "1" (read)
#                       [--changing-between-group-dpsi CHANGING_BETWEEN_GROUP_DPSI] 		--> default = "0.2"
#                       [--non-changing-between-group-dpsi NON_CHANGING_BETWEEN_GROUP_DPSI]
#                       [--changing-between-group-dpsi-secondary CHANGING_BETWEEN_GROUP_DPSI_SECONDARY]--> default="0.1"
#                       [--non-changing-median-reads-threshold NON_CHANGING_MEDIAN_READS_THRESHOLD]
#                       [--permissive-event-non-changing-threshold PERMISSIVE_EVENT_NON_CHANGING_THRESHOLD]
#                       [--non-changing-pvalue-threshold NON_CHANGING_PVALUE_THRESHOLD] 		--> default = "0.05"
#                       [--non-changing-within-group-IQR NON_CHANGING_WITHIN_GROUP_IQR]			--> default = "0.1"
#                       [--changing-pvalue-threshold CHANGING_PVALUE_THRESHOLD]					--> default = "0.05"
#                       [--probability-changing-threshold PROBABILITY_CHANGING_THRESHOLD]		-->default = "0.95"
#                       [--probability-non-changing-threshold PROBABILITY_NON_CHANGING_THRESHOLD]--> default = "0.95"
#                       -d DIRECTORY [-j NPROC] [--debug] [--lazy-load-zarr]
#                       [--parallel-chunksize PARALLEL_CHUNKSIZE] [-l LOGGER]
#                       [--silent]
#                       files [files ...] # splicegraphs and voila files


mkdir -p ${work_dir}/modulizer/
mkdir -p ${work_dir}/modulizer/${CURRENT_SAMPLE}

srun pixi run voila --license ${lic} modulize \
		--nproc ${SLURM_CPUS_PER_TASK} \
		-d ${work_dir}/modulizer/${CURRENT_SAMPLE}/ \
		--show-all --overwrite  \
		${work_dir}/build_output/sg_${CURRENT_SAMPLE}.zarr \
        ${work_dir}/sg_coverage/${CURRENT_SAMPLE}_vs_all.sgc \
        ${work_dir}/sg_coverage/all_vs_${CURRENT_SAMPLE}.sgc \
        ${work_dir}/deltapsi/voila/voila_${CURRENT_SAMPLE}.dpsicov ## contains psi-cov of Sample-vs-all & all-vs-Sample.

       
#--> @voila.log : ERROR:  "CRITICAL - No Majiq License files were found, exiting" 
#--> SOLUTION : "--license" just after <voila>


#> OUTPUT :https://biociphers.bitbucket.io/majiq-docs/modulizer/output.html
# summary.tsv
# cassette.tsv
# tandem_cassette.tsv
# alt5prime.tsv
# alt3prime.tsv
# alt3and5prime.tsv
# etc..



###
### VIEW 
### 


## instructions for "VIEW": https://biociphers.bitbucket.io/majiq/VOILA_view_server.html

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

