## Run all commands
All: Unpack Summary Index Align Counts Coverage

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
	cat ids | parallel "hisat2 --rna-strandness R --max-intronlen 2500 -x refs/grinch-genome.fa -U reads/{}.fq  | samtools sort > bam/{}.bam"
## Create bam index for all files
	cat ids | parallel samtools index bam/{}.bam

## Generate counts
Counts:
## Generate antisense counts
	featureCounts -s 1 -a refs/grinch-annotations_3.gtf -o counts-anti.txt bam/C*.bam bam/W*.bam
## Generate sense counts
	featureCounts -s 2 -a refs/grinch-annotations_3.gtf -o counts-sense.txt bam/C*.bam bam/W*.bam

## Calculate percent coverage
Coverage:
## Select the genes from the GTF file.
	cat refs/grinch-annotations_2.gff | awk '$$3=="gene" { print $$0 }' > genes.gff
## Compute the coverage of each gene relative to all BAM files
	bedtools coverage -S -a genes.gff -b bam/*.bam > coverage.txt
## Sort and store results
	cat coverage.txt | cut -f 9,13 | tr ";" "\t" | cut -f 1,3 | sort -k2,2rn > tin.txt