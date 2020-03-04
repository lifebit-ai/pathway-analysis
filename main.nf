#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/pathway-analysis
========================================================================================
 lifebit-ai/pathway-analysis Nextflow pipeline to perfrom pathway analysis for RNASeq FASTQ files 
 #### Homepage / Documentation
 https://github.com/lifebit-ai/pathway-analysis
----------------------------------------------------------------------------------------
*/

Channel
  .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
  .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
  .set { reads }

Channel
  .fromPath( params.transcriptome )
  .map { file -> tuple(file.baseName, file) }
  .ifEmpty { exit 1, "Cannot find any transcriptome file : ${params.transcriptome}" }
  .set { transcriptome }

/*--------------------------------------------------
  Index the transcriptome
---------------------------------------------------*/

process index {
  tag "$name"

  input:
  set val(name), file(transcriptome) from transcriptome

  output:
  file "${name}.index" into transcriptome_index

  script:
  """
  kallisto index -i ${name}.index ${transcriptome}
  """
}

/*--------------------------------------------------
  Map reads to the transcriptome
---------------------------------------------------*/

process mapping {
  tag "$name"
  publishDir params.outdir

  input:
  file index from transcriptome_index
  set val(name), file(reads) from reads

  output:
  file "kallisto_${name}" into kallisto_out_dirs
  file "stdout.txt" into kallisto_results

  script:
  single_flag = params.singleEnd ? "--single -l ${params.fragment_len} -s ${params.fragment_sd}" : ''
  """
  mkdir kallisto_${name}
  kallisto quant ${single_flag} -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads} &>stdout.txt
  """
}