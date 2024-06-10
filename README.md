# te_resolver

## Goal:
Create pipeline to find TE within genome

## Dictionary:
TE - transposable elements
ONT - Oxford Nanopore Technologies

## Data:
Reference genome - in `genome` directory
Transcposon sequences (TE) - in `te` directory
Long reads from ONT sequencing - in `ont` directory

## Preparation
* Install local blast `sudo apt install ncbi-blast+`
* Install bedtools `sudo apt install bedtools`
* Install seqtk `sudo apt install seqtk`

## How to run
* Do preparations
* Clone te resolver `git clone https://github.com/pkula/te_resolver.git`
* Create directories: `genome, ont, te` for source files
* copy genome in fasta format to genome directory
* copy te to te directory
* copy ont to ont directory
* run python code - the simplest way `python main.py`
