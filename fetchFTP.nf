#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.genomes = false
params.ftp_threads = 100

// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/BOFFO/fetchFTP.nf <ARGUMENTS>
    
    Required Arguments:
      --genomes             CSV file listing genomes (from https://www.ncbi.nlm.nih.gov/genome/browse)
      --output_folder       Folder to write output files to

    Standalone utility to download a set of genomes and create a new manifest CSV which
    points to those downloaded files instead of the original FTP URI.

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || !params.genomes || !params.output_folder){
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

    // Parse the manifest to get a name and FTP prefix for each genome
    sanitize_manifest(
        file(params.genomes)
    )

    // Reformat the manifest to use the updated file paths
    reformat_manifest(
        sanitize_manifest.out
    )

    // Split up the manifest CSV to process each individual genome
    sanitize_manifest.out.map {
        r -> r.splitCsv(
            header: true
        )
    }.flatten(
    ).branch {
        remote: it["uri"] == ""
        local: true
    }.set {
        split_genome_ch
    }

    // Download the genomes from the NCBI FTP server
    fetchFTP(
        split_genome_ch.remote.map {
            r -> [
                r["GenBank FTP"].split("/")[-1].replaceAll(/"/, ""), // Unique genome ID
                r["#Organism Name"],                                 // Readable name
                r["GenBank FTP"]                                     // FTP folder containing the genome
            ]
        }
    )

}

// #############
// # PROCESSES #
// #############

// Fetch genomes via FTP
process fetchFTP {
    tag "Download NCBI genomes by FTP"
    container 'quay.io/fhcrc-microbiome/wget:latest'
    label 'io_limited'
    errorStrategy "retry"
    maxForks params.ftp_threads
    publishDir params.output_folder, mode: 'copy', overwrite: true

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

// Parse the manifest and sanitize the fields
process sanitize_manifest {
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

    input:
        path "raw.manifest.csv"
    
    output:
        path "manifest.csv", emit: manifest
    
"""
#!/usr/bin/env python3

import pandas as pd
import re

df = pd.read_csv("raw.manifest.csv")

print("Subsetting to three columns")
df = df.reindex(
    columns = [
        "GenBank FTP",
        "#Organism Name",
        "uri"
    ]
)

# Remove rows where the "GenBank FTP" doesn't start with "ftp://"
input_count = df.shape[0]
df = df.loc[
    (df["GenBank FTP"].fillna(
        ""
    ).apply(
        lambda n: str(n).startswith("ftp://")
    )) | (
        df["uri"].fillna("").apply(len) > 0
    )
]
print("%d / %d rows have valid FTP or file paths" % (input_count, df.shape[0]))

# Force organism names to be alphanumeric
df = df.apply(
    lambda c: c.apply(lambda n: re.sub('[^0-9a-zA-Z .]+', '_', n)) if c.name == "#Organism Name" else c
)

df.to_csv("manifest.csv", index=None, sep=",")
"""
}

// Format the manifest using the output directory for all FTP files
process reformat_manifest {
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"
    publishDir params.output_folder, mode: 'symlink'

    input:
        path "raw.manifest.csv"
    
    output:
        path "manifest.csv", emit: manifest
    
"""
#!/usr/bin/env python3

import os
import pandas as pd

output_folder = "${params.output_folder}"

# Function to get the 'uri' for each genome
def format_uri(r):

    # Check if the 'uri' is missing
    if r["uri"] is None or isinstance(r["uri"], float) or len(r["uri"]) == 0:

        # Get the UUID inferred from the GenBank FTP path
        genome_uuid = r["GenBank FTP"].split("/")[-1].replace('"', '')

        return os.path.join(
            output_folder,
            f"{genome_uuid}.fasta.gz"
        )

    # The 'uri' is _not_ missing
    else:

        # Return the pre-specified uri
        return r["uri"]

df = pd.read_csv("raw.manifest.csv")

# Add the ultimate file path to the 'uri' column
df = df.assign(
    uri = df.apply(format_uri, axis=1)
)

df.to_csv("manifest.csv", index=None, sep=",")
"""
}
