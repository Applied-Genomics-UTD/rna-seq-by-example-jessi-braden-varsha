## Run all commands
All: Unpack Summary

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