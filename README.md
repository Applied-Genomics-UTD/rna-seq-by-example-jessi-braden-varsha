## Change git repository in RNA-Seq script before running
## Install all necessary packages and environments
## Install bedGraphToBigWig
conda install -y ucsc-bedgraphtobigwig
## Create salmon environment
conda create -y -n salmon
## Install salmon and kallisto
conda install salmon kallisto parallel
## Create a new environment for statistics.
conda create --name stats python=3.8 -y
## Install the statistical packages for rna-seq.
URL=http://data.biostarhandbook.com/books/rnaseq/code/rnaseq-conda.txt
curl $URL | xargs conda install -y
## Run RNA-Seq script
bash RNA-Seq.sh