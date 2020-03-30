# Bacterial Operon Finder for Functional Organization

Tool to find search bacterial genomes for candidate operons

## Goal

Search through a set of bacterial genomes for a set of genes, identify
the location of each gene in each genome, and return a set of tables
which nicely organizes this information for the user.

At the end of the day, a biologist should be able to find which genomes
contain an operon of interest, including the organization of those genes
as well as their protein sequences.

## Approach

Using a Nextflow framework, BOFFO will:

- Align each query protein against each query genome using BLAST
- Collect information for each genome
- Return a set of tables with the level of similarity and organization of each gene in each genome

The input data will be:

- A multi-FASTA with the amino acid sequence of all genes in the operon
- A CSV table listing the genomes from NCBI to query

The CSV table of genomes to query can be obtained from the [NCBI Genome Portal](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/).
This portal allows you to filter to your organisms of interest and
download a summary of those genomes as a CSV table. That is the file
which can be used as an input to BOFFO.

## Usage

```
Usage:

nextflow run FredHutch/BOFFO <ARGUMENTS>

Required Arguments:
  --genomes             CSV file listing genomes (from https://www.ncbi.nlm.nih.gov/genome/browse)
  --operon              Amino acid sequences to search for, in multi-FASTA format
  --output_folder       Folder to write output files to
  --output_prefix       Prefix to use for output file names

Optional Arguments:
  --min_identity        Percent identity threshold used for alignment (default: 90)
  --min_coverage        Percent coverage threshold used for alignment (default: 50)
```
