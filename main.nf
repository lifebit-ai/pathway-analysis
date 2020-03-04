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