# pathway-analysis

Nextflow pipeline to perfrom pathway analysis for RNA-Seq FASTQ files

## Example usage:
```bash
nextflow run main.nf --reads "testdata/*.fastq" --singleEnd --annotation testdata/experiment.csv --transcriptome testdata/transcriptome.fa  --condition treatment
```

