# pathway-analysis

Nextflow pipeline to perfrom pathway analysis for RNA-Seq FASTQ files

## Example usage:
```bash
nextflow run main.nf \
  --experiment testdata/experiment.csv \
  --transcriptome https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/transcriptome/transcriptome.fa \
  --condition treatment
```