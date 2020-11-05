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
  --operon_list         Flag to input multiple sequences per gene -- see description below
  --output_folder       Folder to write output files to
  --output_prefix       Prefix to use for output file names

Optional Arguments:
  --min_identity        Percent identity threshold used for alignment (default: 90)
  --min_coverage        Percent coverage threshold used for alignment (default: 50)
  --max_operon_gap      Maximum gap between genes in the same 'operon' (only used for the 'operon_context' output column) (default: 10000)
  --batchsize           Number of samples to join in each batch (default: 100)
  --max_evalue          Maximum E-value threshold used to filter initial alignments (default: 0.001)
  --annotations         If specified, annotate the regions of all genomes which contain operons
  --annotation_window   The additional area on either side of the operon to annotate (in bp) (default: 10000)

Operon List Input:

Instead of providing a single representative sequence per gene, the --operon_list flag
can be used to point to a comma-delimited list of FASTA files, each containing multiple
sequences from the gene of interest. If this flag is selected, then the conserved positions
within each gene are first identified using PSI-blast to create a position-specific
scoring matrix (PSSM), and then that PSSM is used to query the input genomes. The format
used to link gene names to multi-FASTA files is as follows:

    --operon_list lacZ=test_data/lacZ.fasta,lacY=test_data/lacY.fasta,lacA=test_data/lacA.fasta

    In the example above, the name of each gene is linked to a multi-FASTA with '=', and
    multiple genes are delimited with ','

When this type of analysis is performed, the --operon flag cannot be used, and additional
outputs will be generated including the PSSM generated for each gene.


Annotations:

Note that including the --annotations flag will result in a larger amount of compute being
executed, including the annotation of input genomes using Prokka and the comparison of those 
genomes using clinker to form an interactive visualization (saved in the html/ output folder).

As an alternative approach, we suggest using the entrypoint 'FredHutch/BOFFO/annotate.nf',
which allows the user to annotate just a subset of the results from a set of BOFFO alignments.


Citations:

Seemann T. Prokka: rapid prokaryotic genome annotation
Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063

Gilchrist, C.L.M, 2020. clinker: Easy gene cluster comparison figure generator.

```
