## Make work directory
mkdir work

## Change directory to work
cd work/

## Download golden snidget genome
wget -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.genome.tar.gz

## Unpack genome
tar xzvf golden.genome.tar.gz

## Activate conda environment
ml load anaconda3

## Activate biostars
source activate biostars

## Evaluate the FASTA files
seqkit stats refs/*.fa

