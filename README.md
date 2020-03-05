# pathway-analysis

Nextflow pipeline to perfrom pathway analysis for RNA-Seq FASTQ files

## Example usage:
```bash
# Get testdata
cd testdata
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/transcriptome/transcriptome.fa && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493366.fastq && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493367.fastq && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493368.fastq && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493369.fastq && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493370.fastq && \
wget https://github.com/lifebit-ai/kallisto-nf/raw/master/tutorial/reads/SRR493371.fastq

nextflow run main.nf --reads "testdata/*.fastq" --singleEnd --annotation testdata/experiment.csv --transcriptome testdata/transcriptome.fa  --condition treatment
```