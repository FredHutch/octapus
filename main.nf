#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.output_prefix = false
params.genomes = false
params.operon = false
params.min_identity = 90
params.min_coverage = 50

// Function which prints help message text
def helpMessage() {
    log.info"""
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

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || !params.genomes || !params.output_folder || !params.output_prefix || !params.operon){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // Make sure that the input files exist
    if (file(params.genomes).isEmpty()){
        log.info"""
        The specified --genomes file cannot be found!
        """.stripIndent()
        exit 1
    }
    if (file(params.operon).isEmpty()){
        log.info"""
        The specified --operon file cannot be found!
        """.stripIndent()
        exit 1
    }

    // Parse the manifest to get a name and FTP prefix for each genome
    Channel.from(
        file(params.genomes)
    ).splitCsv(
        header: true
    ).map {
        r -> [
            r["GenBank FTP"].split("/")[-1].replaceAll(/"/, ""), // Unique genome ID
            r["#Organism Name"],                                 // Readable name
            r["GenBank FTP"]                                     // FTP folder containing the genome
        ]
    }.filter {
        it[1].length() > 1
    }.set {
        split_genome_ch
    }

    // Download the genomes from the NCBI FTP server
    fetchFTP(
        split_genome_ch
    )

}

// #############
// # PROCESSES #
// #############

// Fetch genomes via FTP
process fetchFTP {
    tag "Download NCBI genomes by FTP"
    container 'quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), val(ftp_prefix)
    
    output:
        file "${uuid}.fasta.gz"
    
"""
#!/bin/bash
set -e


echo "Downloading ${uuid} from ${ftp_prefix}"

wget --quiet -O ${uuid}.fasta.gz ${ftp_prefix}/${uuid}_genomic.fna.gz

# Make sure the file is gzip compressed
(gzip -t ${uuid}.fasta.gz && echo "${uuid}.fasta.gz is in gzip format") || ( echo "${uuid}.fasta.gz is NOT in gzip format" && exit 1 )

"""
}