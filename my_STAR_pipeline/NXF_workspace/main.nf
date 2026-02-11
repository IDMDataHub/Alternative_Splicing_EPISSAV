#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// ~~~ Subject ~~
// Align sequences of human samples for EPIS-SAV project,
// Count gene-expression and splicing events of neuromuscular diseases

// Author: Maria Kondili
// Date : 08 Augoust 2025, Final Version: 25 Sept 2025



workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}


//params.samplelist   = "/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/my_STAR_pipeline/Samplesheets/all_samples_to_run.tsv" 
// if you want to give in sbatch command, add "null" so the variable is created:
params.samplelist   = null 
params.data_dir     = "/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/raw_data/merged_fastq/"
params.outdir       = "/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/my_STAR_pipeline/Outputs" // also given in nextflow.config !

//
// References
// 
params.starIndex = "/shared/bank/homo_sapiens/GRCh38/star-2.7.5a"  // chr-numb without "chr"
//params.ref_Gtf      = "/shared/home/mkondili/genomes/Annotation/hg38/gencode.v39.annotation.gtf" // chr-numb with "chr" !!!
params.gtf       = "/shared/bank/homo_sapiens/GRCh38.p14/Ensembl_110/gtf/Homo_sapiens.GRCh38.110.gtf" // chr-numb without "chr" ~ in accordance with star-index,and fasta
params.fasta     = "/shared/bank/homo_sapiens/GRCh38.p14/Ensembl_110/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" // chr-numb without "chr"
params.star_path = "/shared/ifbstor1/software/miniconda/envs/star-2.7.5a/bin/STAR"
params.rsem_ref  =  "/shared/projects/splicing_events_detection_rare_diseases/EPIS_SAV/my_STAR_pipeline/NXF_workspace/RSEM_Ref/hg38_transcripts_ref" // created with: prepare_RSEM_ref.sh
//params.cpus = 20 --> global, for all processes. But, inside each process I defined different cpus

//
// create input fastq as variables
//

inputChannel = Channel.fromPath("${params.samplelist}")
    .splitCsv(header:false, sep:'\t') 
    .map { cols ->
        def sample = cols[0].trim()
        def r1 = file(cols[1].trim())
        def r2 = file(cols[2].trim())
        tuple(sample, r1, r2)
    }
	



//
// Define Processes
//


process FastQC_step1 {

  cpus 2
  memory "6G"
  module "fastqc/0.12.1"

  //copy the outputs in an external folder of nextflow, where I want
  publishDir "${params.outdir}/FastQC_before/", mode: 'copy'

  input: 
  tuple val(s), path(r1), path(r2) // tuple = array of 3 values

  output:
  tuple val(s),path("*_fastqc.html")


  shell:
  """
  echo "~~~> FastQC-step1 for !{s} is running now.."
  fastqc  -t ${task.cpus} !{r1} !{r2}

  """
  // ATTENTION: variables from nextflow will be used with "!{}" inside bash, while bash/cluster variables with "${}".
}



process Cut_Adapters {

    cpus 10
    memory "15G"
    module "cutadapt/4.0"

    publishDir "${params.outdir}/CutAdapted", mode: 'copy'

    input:
    tuple val(s), path(r1), path(r2)

    output:
    tuple val(s), path("${s}_R1_cutadapt.fastq.gz"), path("${s}_R2_cutadapt.fastq.gz")

    shell:
	"""

	module li ;
	echo "~~~> CutAdapter will be launched for !{s}"
	cutadapt \
        --cores ${task.cpus} \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minimum-length 100 \
        --error-rate 0.1 \
        --report minimal --discard-trimmed \
        --output !{s}_R1_cutadapt.fastq.gz \
        -p       !{s}_R2_cutadapt.fastq.gz \
        !{r1} !{r2}
    """

    // https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm 
    // Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    // Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
}

process FastQC_step2 {

  cpus 2
  memory "6G"
  module "fastqc/0.12.1"

  publishDir "${params.outdir}/FastQC_After/", mode: 'copy'

  input: 
  tuple val(s), path(r1), path(r2)

  output:
  tuple val(s),path("*_fastqc.html")


  shell:
  """
  echo "~~~> FastQC-step2 for !{s} is running now.."
  fastqc  -t ${task.cpus} !{r1} !{r2}

  """

}


process STAR_Alignment {

      cpus 20
      memory "60G"
      module "star/2.7.5a:perl/5.26.2"

      //copy the outputs in an external folder of nextflow
      publishDir "${params.outdir}/Aligned", mode: 'copy'

      input:
      tuple val(sample),path(f1),path(f2)

      output:
      tuple val(sample),path("${sample}_Aligned.toTranscriptome.out.bam"),path("${sample}_Aligned.sortedByCoord.out.bam")

      shell:
      """
      module li ;
      echo "~~> STAR alignment will be launched now for:"
      echo "SampleID = !{sample}";

      STAR --runMode  alignReads \
            --quantMode TranscriptomeSAM \
            --outSAMattributes NH HI AS NM XS \
            --outSAMtype BAM SortedByCoordinate \
            --quantTranscriptomeBAMcompression 10 \
            --runThreadN ${task.cpus} \
            --genomeDir !{params.starIndex} \
            --sjdbGTFfile !{params.gtf} \
            --readFilesIn !{f1} !{f2}  \
            --outFileNamePrefix !{sample}_  \
            --outReadsUnmapped Fastx  \
            --readFilesCommand zcat
            ##--sjdbFileChrStartEnd   # if junction-regions known to count on
            ## zcat --> use unzipped fastq,when .gz

       """
}

/* ABOUT  SJ-DB-FILE-CHR_START_END: list of annotated junctions

STAR can also utilize annotations formatted as a list of splice junctions coordinates in a text file:
--sjdbFileChrStartEnd /path/to/sjdbFile.txt. This file should contains 4 columns separated
by tabs:
Chr \tab Start \tab End \tab Strand=+/-/.
Here Start and End are first and last bases of the introns (1-based chromosome coordinates). This
file can be used in addition to the --sjdbGTFfile, in which case STAR will extract junctions from
both files


/* ABOUT SAM-Attibutes RSEM expects: 
RSEM parses optional SAM tags to understand multi-mapping and scoring. Common ones:
NH:i:1
→ Number of reported alignments for this read (1 = uniquely mapped).
Critical for RSEM to know how to fractionally assign multi-mappers.
HI:i:1
→ Hit index among multiple alignments for the same read.
(If a read maps 3 times, you’d see HI:i:1, HI:i:2, HI:i:3 across them.)
AS:i:0
→ Alignment score for this alignment (from the aligner).
Used to evaluate quality of mapping.
nM:i:0 (or NM:i:0)
→ Edit distance: number of mismatches between read and reference.
XS:A:+
→ Strand of the RNA-seq alignment (needed for spliced reads in TopHat/STAR).
*/ 



process MarkDuplicates {

	cpus 4
	memory "30GB"
	module "picard/2.22.0"

	publishDir "${params.outdir}/MarkDuplicates/", mode: 'copy'

	input:
		tuple val(sample),path(bam2transcr),path(bam2coord)
	
	output:
		tuple val(sample),path("${sample}_Aligned.sortedByCoord.mrkDups.bam")

	shell:

		"""
		module li;
		echo "~~~> Marking & removing seq. Duplicates from the bam-files now...";
		
		#no need of "java -jar picard.jar" in IFB, just "picard"
		picard -Xmx30g MarkDuplicates \
						REMOVE_SEQUENCING_DUPLICATES=true \
						ASSUME_SORTED=true \
						I=!{bam2coord} \
						O=!{sample}_Aligned.sortedByCoord.mrkDups.bam \
		                M=!{sample}.MarkDuplicates.metrics.txt


		"""
}



process Index_Bam { 
  cpus 16
  memory "32G"
  module "samtools/1.21"

  
  publishDir "${params.outdir}/Aligned", mode: 'copy'

  input:
    tuple val(sample),path(bam_mrkDup) 

  output:
    tuple val(sample),path("${sample}_Aligned.sortedByCoord.out.bam.bai")

  shell:
  	""" 
  	echo "~~> Running Indexing..."
  	samtools index --bai -@ 16 !{bam2coord}

  	"""
}


/*
process RSEM_prepare_ref { --> run in sbatch , only once
	cpus 10
	module "samtools/1.13:rsem/1.3.2:star/2.7.5a"
	input:
		path(star_path)
		path(fasta)
        	path(gtf) 
        
    	output:
    		path("star_hg38_ref.*")
	shell:
		"""
		rsem-prepare-reference --gtf !{gtf} \
		--num-threads $task.cpus \
		!{fasta}  ./RSEM/hg38_transcript_ref
		"""
}
*/


process RSEM_Quantification {

	cpus 16
	memory "32G"
	module "samtools/1.21:rsem/1.3.2:star/2.7.5a"

	
	publishDir "${params.outdir}/RSEM_counts", mode: 'copy'

	input:
	   tuple val(sample),path(bam2transcr),path(bam2coord)

	output:
	   tuple val(sample),
           path("${sample}.isoforms.results"),
           path("${sample}.genes.results")

	shell:
 	""" 
    echo "~~> Running RSEM quantification for sample !{sample}"

	rsem-calculate-expression \
        --strandedness reverse \
        --num-threads $task.cpus \
        --no-bam-output \
        --alignments \
        --paired-end  \
        !{bam2transcr} \
        !{params.rsem_ref} \
        !{sample}
        """
}




//
// Call Workflow modules 
//

workflow {

	inputChannel.view{ sample, r1, r2 -> "~~>Sample ${sample},\nR1 :${r1},\nR2: ${r2}" }

	qcBeforeChannel = FastQC_step1(inputChannel)

	cutadaptChannel = Cut_Adapters(inputChannel) 

	qcAfterChannel  = FastQC_step2(cutadaptChannel)

	alignedChannel  = STAR_Alignment(cutadaptChannel)

	markDupChannel  = MarkDuplicates(alignedChannel)

	indexedChannel  = Index_Bam(markDupChannel) // to use mrkDup.bam for FRASER, they should be indexed

	rsemChannel     = RSEM_Quantification(alignedChannel) // | view

	// RSEM only for Aligned.toTranscriptome
}



