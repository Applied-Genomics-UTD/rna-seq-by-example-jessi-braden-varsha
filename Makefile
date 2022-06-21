## Run all commands
All: Unpack Summary Index Align

## Download grouchy grinch data
grinch.tar.gz:
	wget -nc http://data.biostarhandbook.com/rnaseq/data/grinch.tar.gz 

## Unpack the data
Unpack: grinch.tar.gz
	tar zxvf grinch.tar.gz

## Make data and genome statistics
Summary:
	seqkit stats reads/*fq > Reads_Stats
	seqkit stats refs/grinch-genome.fa > Genome_Stats

## Create ids file
ids:
	parallel echo {1}{2} ::: Cranky Wicked ::: 1 2 3 > ids

## Index genome
Index:
	hisat2-build refs/grinch-genome.fa refs/grinch-genome.fa

## Run aligner in parallel
Align: ids
## Make bam directory
	mkdir -p bam
## Run hisat2 for all samples
	cat ids | parallel "hisat2 --max-intronlen 2500 -x refs/grinch-genome.fa -U reads/{}.fq  | samtools sort > bam/{}.bam"
## Create bam index for all files
	cat ids | parallel samtools index bam/{}.bam