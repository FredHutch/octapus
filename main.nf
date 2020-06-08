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
params.max_operon_gap = 10000
params.batchsize = 100

// Import modules
include {
    collectResults as collectResultsRound1;
    collectResults as collectResultsRound2;
    collectFinalResults;
    summaryPDF;
} from './modules/modules' params(
    output_prefix: params.output_prefix,
    output_folder: params.output_folder,
)

// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__plotting = "quay.io/fhcrc-microbiome/boffo-plotting:latest"


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
      --max_operon_gap      Maximum gap between genes in the same 'operon' (only used for the 'operon_context' output column) (default: 10000)
      --batchsize           Number of samples to join in each batch (default: 100)

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
    sanitize_manifest(
        file(params.genomes)
    )
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

    // Align the operon against each genome
    runBLAST(
        fetchFTP.out.mix(
            split_genome_ch.local.map {
                r -> [
                    r["uri"].split("/")[-1],
                    r["#Organism Name"],
                    file(r["uri"])
                ]
            }
        ),
        operon_fasta
    )

    // Parse each individual alignment
    parseAlignments(
        runBLAST.out
    )

    // Join the parsed alignments with the FASTA for each genome
    // For each alignment, extract the sequence of the aligned region
    extractAlignments(
        fetchFTP.out.join(
            parseAlignments.out.map {
                r -> [r.name.replaceAll(/.csv.gz/, ""), r]
            }
        )
    )

    // Collect results in rounds
    collectResultsRound1(
        extractAlignments.out.collate(params.batchsize)
    )
    collectResultsRound2(
        collectResultsRound1.out.collate(params.batchsize)
    )

    // Make a single output table
    collectFinalResults(
        collectResultsRound2.out.collect()
    )

    // Make the FASTA of aligned sequences
    formatFASTA(
        collectFinalResults.out
    )

    // Make the summary PDF
    summaryPDF(
        collectFinalResults.out
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

// Parse each individual alignment file
process parseAlignments {
    tag "Identify operons"
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), path(aln_gz)
    
    output:
        file("${uuid}.csv.gz") optional true
    
"""
#!/usr/bin/env python3
import pandas as pd
import json
import gzip

# Function to figure out the "genome_context" for each gene
def genome_context(df):

    # Return the list of genes detected in the genome, relative to each alignment
    output = []

    # Make the list for the forward strand
    fwd_str = " :: ".join([
        "%s (%s)" % (i["gene_name"], i["strand"])
        for _, i in df.iterrows()
    ])
    # and the reverse strand
    rev_str = " :: ".join([
        "%s (%s)" % (i["gene_name"], "-" if i["strand"] == "+" else "+")
        for _, i in df.iterrows()
    ][::-1])

    # For each gene, use the appropriate list
    for ix, r in df.iterrows():
        if r["strand"] == '+':
            output.append(fwd_str)
        # For genes on the - strand, switch the strand reported for the other genes as well as the order
        else:
            output.append(rev_str)

    return pd.Series(output, index=df.index)

# Function to figure out the "operon_context" for each gene
def operon_context(df, maximum_gap):

    # Walk through the genome and find each of the "operons"
    operons = []
    last_contig = None
    last_pos = None
    for ix, r in df.iterrows():

        # Start of a new contig, and therefore operon
        if last_contig is None or r["contig_name"] != last_contig:
            operons.append([ix])

        # Add to an existing operon
        elif (min(r["contig_start"], r["contig_end"]) - last_pos) < ${params.max_operon_gap}:
            operons[-1].append(ix)
        
        # Outside the maximum gap, so start a new operon
        else:
            operons.append([ix])

        # Set the temp values for the next iteration
        last_contig = r["contig_name"]
        last_pos = max(r["contig_start"], r["contig_end"])

    # Return the list of genes detected in the genome, relative to each alignment, within a given distance
    output = []

    # Walk through each operon
    for operon in operons:
        # Format a string for the positive strand
        fwd_str = " :: ".join([
            "%s (%s)" % (i["gene_name"], i["strand"])
            for _, i in df.reindex(index=operon).iterrows()
        ])
        # and the reverse strand
        rev_str = " :: ".join([
            "%s (%s)" % (i["gene_name"], "-" if i["strand"] == "+" else "+")
            for _, i in df.reindex(index=operon).iterrows()
        ][::-1])

        for ix in operon:
            if df.loc[ix, "strand"] == "+":
                output.append(fwd_str)
            else:
                output.append(rev_str)
            
    return pd.Series(output, index=df.index)

print("Processing alignments for ${genome_name.replaceAll(/"/, "")} against ${uuid}")

# Read in the alignment table
df = pd.read_csv(
    "${aln_gz}", 
    compression = "gzip", 
    sep = "\\t",
    header = None,
    names = [
        "gene_name",
        "contig_name",
        "pct_iden",
        "alignment_length",
        "mismatch",
        "gapopen",
        "gene_start",
        "gene_end",
        "gene_len",
        "contig_start",
        "contig_end",
        "contig_len"
    ]
)

# Calculate coverage and filter by coverage and identity
df = df.assign(
    gene_cov = 100 * ((df["gene_end"] - df["gene_start"]) + 1).abs() / df["gene_len"]
).query(
    "gene_cov >= ${params.min_coverage}"
).query(
    "pct_iden >= ${params.min_identity}"
)

# Check to see if there are any alignments passing the filter
if df.shape[0] == 0:
    print("Zero alignments found passing the filter -- skipping")

else:


    # Assign the strand for each gene
    df = df.assign(
        strand = df.apply(lambda r: "+" if r["contig_start"] < r["contig_end"] else "-", axis=1)
    )

    # Sort by contig length (descending) and alignment position (ascending)
    df.sort_values(
        by = ["contig_len", "contig_name", "contig_start"],
        ascending = [False, True, True],
        inplace = True
    )

    # Figure out the "genome_context" and "operon_context" for each hit
    df = df.assign(
        genome_context = genome_context(df),
        operon_context = operon_context(df, ${params.max_operon_gap}),
    )

    # Add the genome ID, name, and operon size
    df = df.assign(
        genome_id = "${uuid}",
        genome_name = "${genome_name.replaceAll(/"/, "")}",
        operon_size = df["operon_context"].apply(lambda n: len(n.split(" :: ")))
    )

    # Write out as a formatted CSV
    print("Writing out to CSV")
    df.reindex(
        columns = [
            "genome_name",
            "genome_id",
            "contig_name",
            "contig_len",
            "gene_name",
            "gene_len",
            "contig_start",
            "contig_end",
            "genome_context",
            "operon_context",
            "operon_size",
            "pct_iden",
            "gene_cov",
            "mismatch",
            "gapopen",
            "strand",
            "alignment_length",
            "gene_start",
            "gene_end",
        ]
    ).to_csv(
        "${uuid}.csv.gz",
        index = None,
        compression = "gzip",
        sep = ","
    )

print("Done")

"""
}


// Extract the sequence of the aligned regions
process extractAlignments {
    container 'quay.io/fhcrc-microbiome/biopython-pandas:latest'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), path(fasta_gz), path(aln_gz)
    
    output:
        path "${uuid}.withSeq.csv.gz"
    
"""
#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import gzip
import pandas as pd

# Read in all of the sequences for this genome
genome_seqs = dict()
for header, seq in SimpleFastaParser(gzip.open(
    "${fasta_gz}",
    "rt"
)):
    # Parse the contig name from the header
    n = header.split(" ")[0]
    print("Read in %s - %d bps" % (n, len(seq)))
    genome_seqs[n] = seq

# Read in the table
df = pd.read_csv("${aln_gz}")

# Define a function to parse the aligned region for a given row
def extract_aligned_sequence(r):

    # Get the information for this row of the table
    contig_name = r["contig_name"]
    start_pos = r["contig_start"]
    end_pos = r["contig_end"]

    # Make sure that this contig is in the FASTA
    assert contig_name in genome_seqs, "Couldn't find %s in FASTA" % contig_name

    # If the alignment is +, return this region
    if start_pos < end_pos:
        return genome_seqs[
            contig_name
        ][start_pos - 1: end_pos]

    # Otherwise return the reverse complement
    else:
        return str(Seq(
            genome_seqs[
                contig_name
            ][end_pos - 1: start_pos]
        ).reverse_complement())

# Define a function to translate the aligned region for a given row
def translate_aligned_sequence(r):

    return str(Seq(
        r["aligned_sequence"]
    ).translate(table=11))

# Add the nucleotide sequences to the table
print("Assigning sequences for each alignment")
df = df.assign(
    aligned_sequence = df.apply(extract_aligned_sequence, axis=1)
)
print("Done")

# Translate the aligned sequence
print("Translating sequences for each alignment")
df = df.assign(
    translated_sequence = df.apply(translate_aligned_sequence, axis=1)
)
print("Done")

# Write out to a file
fpo = "${uuid}.withSeq.csv.gz"
print("Writing out to %s" % fpo)
df.to_csv(fpo, index=None)
print("Done")
"""
}


// Format FASTA files with aligned sequences
process formatFASTA {
    container 'quay.io/fhcrc-microbiome/biopython-pandas:latest'
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true

    input:
        file results_csv_gz
    
    output:
        path "*/*.gz"
    
"""
#!/usr/bin/env python3

import gzip
import os
import pandas as pd

# Read in all of the results for this dataset
df = pd.read_csv("${results_csv_gz}")

# For each operon structure, make a folder
for operon_structure, operon_df in df.groupby("operon_context"):

    # Format the folder name
    operon_folder = operon_structure.replace(
        "::", "_"
    ).replace(
        "(+)", "FWD"
    ).replace(
        "(-)", "REV"
    ).replace(
        " ", "_"
    ).replace(
        "__", "_"
    ).replace(
        "__", "_"
    )

    print("Writing out %d sequences for %s" % (operon_df.shape[0], operon_structure))
    os.mkdir(operon_folder)

    # For each gene, make a FASTA file with the nucleotide sequence
    for gene_name, gene_df in df.groupby("gene_name"):
        with gzip.open("%s/%s.alignments.fasta.gz" % (operon_folder, gene_name), "wt") as fo:
            fo.write("\\n".join([
                ">%s::%s::%s::%s-%s\\n%s" % (
                    r["genome_id"], 
                    r["contig_name"],
                    r["genome_name"],
                    r["contig_start"],
                    r["contig_end"],
                    r["aligned_sequence"],
                )
                for _, r in gene_df.iterrows()
            ]))

        # Also write out the translated sequence
        with gzip.open("%s/%s.alignments.fastp.gz" % (operon_folder, gene_name), "wt") as fo:
            fo.write("\\n".join([
                ">%s::%s::%s::%s-%s\\n%s" % (
                    r["genome_id"], 
                    r["contig_name"],
                    r["genome_name"],
                    r["contig_start"],
                    r["contig_end"],
                    r["translated_sequence"],
                )
                for _, r in gene_df.iterrows()
            ]))
"""
}
