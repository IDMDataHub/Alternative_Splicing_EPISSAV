#!/usr/bin/env bash

cd /shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/rMATS/input_bam/

#cat sample_*.csv | grep -v "C002I3K" | paste -sd, - > sample_K_vs_all.csv #ok
cat sample_*.csv | grep -v "C002I44" | paste -sd, - > sample_44_vs_all.csv
cat sample_*.csv | grep -v "C002I45" | paste -sd, - > sample_45_vs_all.csv
cat sample_*.csv | grep -v "C002I47" | paste -sd, - > sample_47_vs_all.csv
cat sample_*.csv | grep -v "C002I4F" | paste -sd, - > sample_4F_vs_all.csv
cat sample_*.csv | grep -v "C002I4G" | paste -sd, - > sample_4G_vs_all.csv
cat sample_*.csv | grep -v "C002I3L" | paste -sd, - > sample_L_vs_all.csv
cat sample_*.csv | grep -v "C002I3M" | paste -sd, - > sample_M_vs_all.csv
cat sample_*.csv | grep -v "C002I3N" | paste -sd, - > sample_N_vs_all.csv
cat sample_*.csv | grep -v "C002I3P" | paste -sd, - > sample_P_vs_all.csv
cat sample_*.csv | grep -v "C002I3Q" | paste -sd, - > sample_Q_vs_all.csv
cat sample_*.csv | grep -v "C002I3R" | paste -sd, - > sample_R_vs_all.csv
cat sample_*.csv | grep -v "C002I3S" | paste -sd, - > sample_S_vs_all.csv
cat sample_*.csv | grep -v "C002I3U" | paste -sd, - > sample_U_vs_all.csv
cat sample_*.csv | grep -v "C002I3W" | paste -sd, - > sample_W_vs_all.csv
cat sample_*.csv | grep -v "C002I3Y" | paste -sd, - > sample_Y_vs_all.csv
cat sample_*.csv | grep -v "C002I3Z" | paste -sd, - > sample_Z_vs_all.csv
