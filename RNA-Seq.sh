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

## Make an index file to download genome into IGV
samtools faidx refs/genome.fa

## Give instructions for next part of analysis
echo "Download features.gff, genome.fa, and genome.fa.fai"
echo "Upload genome.fa and genome.fa.fai to <Genome> in IGV at the same time"
echo "Upload features.gff to <Tracks> in IGV"

## Download golden snidget reads
wget -nc http://data.biostarhandbook.com/books/rnaseq/data/golden.reads.tar.gz

## Unpack reads
tar zxvf golden.reads.tar.gz

## Run statistics on reads
seqkit stats reads/*.fq > read_stats

## Generate root names with parallel and save to file "ids"
parallel -j 1 echo {1}_{2} ::: BORED EXCITED ::: 1 2 3 > ids

## Set and export reference genome
export IDX=refs/genome.fa

## Build the genome index
hisat2-build $IDX $IDX

## Index the reference with samtools
samtools faidx $IDX

## Create BAM folder
mkdir -p bam

## Align FASTQ files to the reference genome
cat ids | parallel "hisat2 -x $IDX -1 reads/{}_R1.fq -2 reads/{}_R2.fq | samtools sort > bam/{}.bam"

## Index each BAM file
cat ids | parallel  "samtools index bam/{}.bam"

## Install bedGraphToBigWig
conda install -y ucsc-bedgraphtobigwig

## Turn each BAM file into bedGraph coverage. The files will have the .bg extension.
cat ids | parallel "bedtools genomecov -ibam  bam/{}.bam -split -bg  > bam/{}.bg"

## Convert each bedGraph coverage into bigWig coverage. The files will have the .bw extension.
cat ids | parallel "bedGraphToBigWig bam/{}.bg  ${IDX}.fai bam/{}.bw"

## Give instructions for next portion of analysis
echo "Download all bigWig files and upload to IGV"

## Run featureCounts for all BAM files
featureCounts -p -a refs/features.gff -o counts.txt bam/BORED_?.bam bam/EXCITED_?.bam

## Show the first five lines of the counts.txt file.
cat counts.txt | head -5

## Create salmon environment
conda create -y -n salmon

## Activate salmon environment
conda activate salmon

## Install salmon and kallisto
conda install salmon kallisto parallel

## Set and export the shortcuts to reference and index.
export REF=refs/transcripts.fa
export IDX1=refs/kallisto.idx
export IDX2=refs/salmon.idx

## Build index with kallisto.
kallisto index -i ${IDX1} ${REF}

## Build the index with salmon.
salmon index -t ${REF} -i ${IDX2}

# Kallisto quantification.
kallisto quant -i ${IDX1} -o BORED_1 reads/BORED_1_R1.fq reads/BORED_1_R2.fq

# Salmon quantification.
salmon quant -i ${IDX2} -l A --validateMappings -1 reads/BORED_1_R1.fq -2 reads/BORED_1_R2.fq  -o out2/BORED_1

## Check counts
cat BORED_1/abundance.tsv | head -5

## Create parent directories
mkdir -p out1 out2

## The name of the index.
export REF=refs/transcripts.fa
export IDX1=refs/kallisto.idx
export IDX2=refs/salmon.idx

## Run Kallisto to classify the reads.
cat ids | parallel kallisto quant -i ${IDX1} -o out1/{} reads/{}_R1.fq reads/{}_R2.fq

## Run Salmon to classify the reads.
cat ids | parallel salmon quant -i ${IDX2} -l A --validateMappings -o out2/{} -1 reads/{}_R1.fq -2 reads/{}_R2.fq

## Create bin directory
mkdir ~/bin

## Download the custom script to combine kallisto or salmon outputs.
curl http://data.biostarhandbook.com/books/rnaseq/code/combine.py > ~/bin/combine

## Make the script executable.
chmod +x ~/bin/combine

## Combine kallisto abundances into single count matrix
cat ids | combine out1 > counts1.txt

## Combine salmon abundances into single count matrix
cat ids | combine out2 > counts2.txt

## Install all required packages
conda install bioconductor-edger bioconductor-deseq2 r-gplots -y