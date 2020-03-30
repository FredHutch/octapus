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

    // Point to the operon file
    operon_fasta = file(params.operon)

    // Make sure that the input files exist
    if (file(params.genomes).isEmpty()){
        log.info"""
        The specified --genomes file cannot be found!
        """.stripIndent()
        exit 1
    }
    if (operon_fasta.isEmpty()){
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

    // Align the operon against each genome
    runBLAST(
        fetchFTP.out,
        operon_fasta
    )

    // Parse each individual alignment
    parseAlignments(
        runBLAST.out
    )

    // Make a single output table
    collectResults(
        parseAlignments.out.toSortedList()
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
        tuple val(uuid), val(genome_name), file("${uuid}.fasta.gz")
    
"""
#!/bin/bash
set -e


echo "Downloading ${uuid} from ${ftp_prefix}"

wget --quiet -O ${uuid}.fasta.gz ${ftp_prefix}/${uuid}_genomic.fna.gz

# Make sure the file is gzip compressed
(gzip -t ${uuid}.fasta.gz && echo "${uuid}.fasta.gz is in gzip format") || ( echo "${uuid}.fasta.gz is NOT in gzip format" && exit 1 )

"""
}

// Align the operon against each individual genome
process runBLAST {
    tag "Align operon"
    container 'quay.io/fhcrc-microbiome/blast@sha256:1db09d0917e52913ed711fcc5eb281c06d0bb632ec8cd5a03610e2c3377e1753'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), path(fasta_gz)
        file operon_fasta
    
    output:
        tuple val(uuid), val(genome_name), file("${uuid}.aln.gz")
    
"""
#!/bin/bash
set -e

echo "Running alignment for ${uuid}"

ls -lahtr

tblastn \
    -query ${operon_fasta} \
    -subject <(gunzip -c ${fasta_gz}) \
    -out ${uuid}.aln \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen" \
    -evalue 0.001

echo "Compressing alignment file"

gzip ${uuid}.aln

echo Done
"""
}

// Parse each individual alignment file
process parseAlignments {
    tag "Identify operons"
    container 'quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), path(aln_gz)
    
    output:
        file "${uuid}.json.gz"
    
"""
#!/usr/bin/env python3
import pandas as pd
import json
import gzip

print("Processing alignments for ${genome_name.replaceAll(/"/, "")} against ${uuid}")

# Read in the alignment table
df = pd.read_csv(
    "${aln_gz}", 
    compression = "gzip", 
    sep = "\\t",
    header = None,
    names = [
        "qaccver",
        "saccver",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "sstart",
        "send",
        "slen"
    ]
)

# Calculate coverage and filter by coverage and identity
df = df.assign(
    qcov = 100 * ((df["qend"] - df["qstart"]) + 1).abs() / df["qlen"]
).query(
    "qcov >= ${params.min_coverage}"
).query(
    "pident >= ${params.min_identity}"
)

# Format a dict with the locations of all genes

output = dict()
output["genome_id"] = "${uuid}"
output["genome_name"] = "${genome_name.replaceAll(/"/, "")}"

for gene_name, gene_df in df.groupby("qaccver"):
    print("Found %s alignments for %s" % (gene_df.shape[0], gene_name))
    output[
        gene_name
    ] = "; ".join([
        "%s: %s - %s" % (r["saccver"], r["sstart"], r["send"])
        for _, r in gene_df.iterrows()
    ])

with gzip.open("${uuid}.json.gz", "wt") as handle:
    json.dump(
        output,
        handle
    )

print("Done")

"""
}

// Parse each individual alignment file
process collectResults {
    tag "Make a single table"
    container 'quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed'
    label 'io_limited'
    publishDir "${params.output_folder}"
    // errorStrategy "retry"

    input:
        file json_gz_list
    
    output:
        file "${params.output_prefix}.csv.gz"
    
"""
#!/usr/bin/env python3
import pandas as pd
import json
import gzip

dat = []

for fp in "${json_gz_list}".split(" "):
    print("Reading in %s" % fp)
    dat.append(
        json.load(
            gzip.open(
                fp,
                "rt"
            )
        )
    )

print("Making a single output table")
df = pd.DataFrame(dat)

print("Writing out to ${params.output_prefix}.csv.gz")
df.to_csv("${params.output_prefix}.csv.gz", index=None, compression="gzip")
print("Done")
"""
}