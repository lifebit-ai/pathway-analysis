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

Channel
  .fromPath( params.annotation )
  .ifEmpty { exit 1, "Cannot find annotation/experiment file : ${params.annotation}" }
  .set { annotation }

Channel
  .fromPath( params.rmarkdown )
  .ifEmpty { exit 1, "Cannot find R Markdown file : ${params.rmarkdown}" }
  .set { rmarkdown }

Channel
  .fromPath( params.hallmark_pathways )
  .ifEmpty { exit 1, "Cannot find Hallmark pathways file : ${params.hallmark_pathways}" }
  .set { hallmark_pathways }
Channel
  .fromPath( params.kegg_pathways )
  .ifEmpty { exit 1, "Cannot find KEGG pathways file : ${params.kegg_pathways}" }
  .set { kegg_pathways }
Channel
  .fromPath( params.mir_pathways )
  .ifEmpty { exit 1, "Cannot find miR targets file : ${params.mir_pathways}" }
  .set { mir_pathways }
Channel
  .fromPath( params.go_pathways )
  .ifEmpty { exit 1, "Cannot find GO annotations file : ${params.go_pathways}" }
  .set { go_pathways }

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
  // TODO: add `transcriptome_index` param & `when` to skip this
  """
  kallisto index -i ${name}.index ${transcriptome}
  """
}

/*--------------------------------------------------
  Map reads to the transcriptome
---------------------------------------------------*/

process mapping {
  tag "$name"
  publishDir params.outdir, mode: 'copy'

  input:
  each file(index) from transcriptome_index
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

/*--------------------------------------------------
  Run differential gene expression analysis 
---------------------------------------------------*/

process gene_expression {
  tag "$annotation"
  publishDir params.outdir, mode: 'copy'

  input:
  file('kallisto/') from kallisto_out_dirs.collect()
  file(annotation) from annotation
  file(rmarkdown) from rmarkdown
  file(hallmark) from hallmark_pathways
  file(kegg) from kegg_pathways
  file(mir) from mir_pathways
  file(go) from go_pathways

  output:
  file("{MultiQC,diffexpr-results.csv}") into results

  script:
  // TODO: pass arg for `kallisto_dir to Rmd
  """
  # copy the rmarkdown into the pwd
  cp $rmarkdown tmp && mv tmp $rmarkdown
  R -e "rmarkdown::render('${rmarkdown}', params = list(annotation='${annotation}',condition='${params.condition}',hallmark='${hallmark}',kegg='${kegg}',mir='${mir}',go='${go}'))"
  mkdir MultiQC && mv DE_with_DEseq2.html MultiQC/multiqc_report.html
  """
}  