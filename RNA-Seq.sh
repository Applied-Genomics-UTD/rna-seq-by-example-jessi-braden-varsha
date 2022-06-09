## Change directory to home
cd

## Change directory to scratch
cd scratch/

## Make git repository
git clone https://github.com/Applied-Genomics-UTD/rna-seq-by-example-jessi-braden-varsha.git

## Change directory to RNA-seq by example
cd rna-seq-by-example-jessi-braden-varsha

## Activate conda environment
ml load anaconda3

## Activate biostars
source activate biostars

## Make work directory
mkdir work

## Change directory to work
cd work/

## Download golden snidget genome
wget -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.genome.tar.gz

## Unpack genome
tar xzvf golden.genome.tar.gz

## Evaluate the FASTA files
seqkit stats refs/*.fa

## What does the genome file look like?
cat refs/genome.fa | head -5

## What is the annotation file?
cat refs/features.gff | head -5

## How many exons (approximate count)?
cat refs/features.gff | grep exon | wc -l

## What do the transcripts look like?
cat refs/transcripts.fa | head -5

## Print out more gene names:
cat refs/transcripts.fa | grep ">" | head
