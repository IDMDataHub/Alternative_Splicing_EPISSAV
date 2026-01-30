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


##--> Input
all_samples="${workDir}/input_bam/all_samples.csv"


outDir="${workDir}/Outputs/one_group_stat_off/output_all_statoff"; 
mkdir -p ${outDir}
tmpDir="${outDir}/tmp_all_statoff/"; 
mkdir -p ${tmpDir}


###---- Run All samples "statOFF" -----

## python ${rmats_dir}/rmats.py \
# rmats.py --b1 ${all_samples} \
# --gtf ${gtf} --readLength 101 \
# --od ${outDir} --tmp ${tmpDir} \
# -t "paired"  --libType fr-firststrand \
# --nthread ${cores} \
# --statoff \
# --task both 


###------ Run STAT step for "1-vs-all" ------

mkdir -p ${outDir}/Stats_out/4F_vs_all

python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
		--new-output-dir ${outDir}/Stats_out/4F_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 0 \
                --group-2-indices 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/4F_vs_all --tmp ${outDir}/Stats_out/tmp_4F_vs_all --task stat

mkdir -p ${outDir}/Stats_out/44_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
		--new-output-dir ${outDir}/Stats_out/44_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 1 \
                --group-2-indices 0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/44_vs_all --tmp ${outDir}/Stats_out/tmp_44_vs_all --task stat

mkdir -p ${outDir}/Stats_out/45_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
		--new-output-dir ${outDir}/Stats_out/45_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 2 \
                --group-2-indices 0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/45_vs_all --tmp ${outDir}/Stats_out/tmp_45_vs_all --task stat


mkdir -p ${outDir}/Stats_out/47_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
		--new-output-dir ${outDir}/Stats_out/47_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 3 \
                --group-2-indices 0,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/47_vs_all --tmp ${outDir}/Stats_out/tmp_47_vs_all --task stat


mkdir -p ${outDir}/Stats_out/4G_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
		--new-output-dir ${outDir}/Stats_out/4G_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 4 \
                --group-2-indices 0,1,2,3,5,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/4G_vs_all --tmp ${outDir}/Stats_out/tmp_4G_vs_all --task stat



mkdir -p ${outDir}/Stats_out/K_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/K_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 5 \
                --group-2-indices 0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/K_vs_all --tmp ${outDir}/Stats_out/tmp_K_vs_all --task stat


mkdir -p ${outDir}/Stats_out/L_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/L_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 6 \
                --group-2-indices 0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/L_vs_all --tmp ${outDir}/Stats_out/tmp_L_vs_all --task stat


mkdir -p ${outDir}/Stats_out/M_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/M_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 7 \
                --group-2-indices 0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/M_vs_all --tmp ${outDir}/Stats_out/tmp_M_vs_all --task stat


mkdir -p ${outDir}/Stats_out/N_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/N_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 8 \
                --group-2-indices 0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/N_vs_all --tmp ${outDir}/Stats_out/tmp_N_vs_all --task stat


mkdir -p ${outDir}/Stats_out/P_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/P_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 9 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/P_vs_all --tmp ${outDir}/Stats_out/tmp_P_vs_all --task stat


mkdir -p ${outDir}/Stats_out/Q_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/Q_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 10 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16

python rmats.py --od ${outDir}/Stats_out/Q_vs_all --tmp ${outDir}/Stats_out/tmp_Q_vs_all --task stat



mkdir -p ${outDir}/Stats_out/R_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/R_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 11 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16
python rmats.py --od ${outDir}/Stats_out/R_vs_all --tmp ${outDir}/Stats_out/tmp_R_vs_all --task stat



mkdir -p ${outDir}/Stats_out/S_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/S_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 12 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16
python rmats.py --od ${outDir}/Stats_out/S_vs_all --tmp ${outDir}/Stats_out/tmp_S_vs_all --task stat


mkdir -p ${outDir}/Stats_out/U_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/U_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 13 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16

python rmats.py --od ${outDir}/Stats_out/U_vs_all --tmp ${outDir}/Stats_out/tmp_U_vs_all --task stat


mkdir -p ${outDir}/Stats_out/W_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/W_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 14 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16
        
python rmats.py --od ${outDir}/Stats_out/W_vs_all --tmp ${outDir}/Stats_out/tmp_W_vs_all --task stat


mkdir -p ${outDir}/Stats_out/Y_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/Y_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 15 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16

python rmats.py --od ${outDir}/Stats_out/Y_vs_all --tmp ${outDir}/Stats_out/tmp_Y_vs_all --task stat


mkdir -p ${outDir}/Stats_out/Z_vs_all
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py \
                --new-output-dir ${outDir}/Stats_out/Z_vs_all \
                --old-output-dir ${outDir} \
                --group-1-indices 16 \
                --group-2-indices 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15

python rmats.py --od ${outDir}/Stats_out/Z_vs_all --tmp ${outDir}/Stats_out/tmp_Z_vs_all --task stat
