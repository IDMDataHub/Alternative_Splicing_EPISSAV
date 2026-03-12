
## Preliminary Analysis


The preliminary analysis consists of quality control of the raw data and alignment of the sequences to a reference genome, here human.

The corresponding code can be found in the "my_STAR_pipeline" folder.\
It is a nextflow pipeline, with:
1. main.nf script containing the processes,
2. nextflow.config containing the configuration relative to the cluster where it will run.
3. the script to execute it on a cluster with slurm-scheduler : "run_nxf_star_rsem_pipe.sbatch"


## Alternative splicing analysis with rMATS


To start analysis with rMATS, two scripts are necessary:

1. run_rmats_statoff_task_both.sh
2. run_rmats_stats_1-vs-all.sbatch 

The first script detects the splicing events for all samples, using the "input_bam/" files. 
The output is directed to "output_all_statoff/" folder.

The second script consists of the statistical analysis of the comparison between samples. \
It uses the output of the first and creates a new folder (Stats_out/X_vs_all/, where X is one sample). \
All samples are compared, each -vs- all others in consecutive commands.


## Downstream Analysis


The Downstream analysis is run with the following scripts: 

1. manip_rMATS_output.Rmd : for merging the outputs of rMATS that are given per splicing-event.

2. Filter_signif_ASE_DiseaseGenes.Rmd : for the filtering of results based on FDR, deltapsi and disease-genes

For the filtering by disease-genes, a folder named "GeneTables_Bonne-Rivier/" is also used, but not published in this repository due to its size.

It consists of multiple excel tables, one for each known disease.

They can be downloaded in this site: [muscle-gene-tables](musclegenetable.fr)


