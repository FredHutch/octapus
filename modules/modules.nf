// Parse the manifest and sanitize the fields
process sanitize_manifest {
    container "${params.container__pandas}"
    label 'io_limited'

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
process collectResults {
    
    container "${params.container__pandas}"
    label 'io_limited'

    input:
        file "input/*"
    
    output:
        file "*.csv.gz"
    
"""
#!/usr/bin/env python3
import gzip
import json
import os
import pandas as pd
from shutil import copyfile

# Parse the list of files to join
csv_gz_list = [
    os.path.join("input", fp)
    for fp in os.listdir("input")
]
assert len(csv_gz_list) > 0

# Make sure all inputs are .csv.gz
for fp in csv_gz_list:
    assert fp.endswith(".csv.gz"), fp

# The output table will be named for one of the inputs
output_fp = csv_gz_list[0].replace("input/", "")

# If there is only one file, just copy it to the output
if len(csv_gz_list) == 1:
    copyfile(
        csv_gz_list[0],
        output_fp
    )
else:

    print("Making a single output table")
    df = pd.concat([
        pd.read_csv(fp)
        for fp in csv_gz_list
    ], sort=True)

    print("Writing out to %s" % output_fp)
    df.reindex(
        columns = [
            "operon_context",
            "operon_size",
            "operon_ix",
            "genome_context",
            "genome_id",
            "genome_name",
            "contig_name",
            "contig_start",
            "contig_end",
            "contig_len",
            "strand",
            "alignment_length",
            "gene_name",
            "gene_start",
            "gene_end",
            "gene_len",
            "gene_cov",
            "gapopen",
            "mismatch",
            "pct_iden",
            "aligned_sequence",
            "translated_sequence",
        ]
    ).to_csv(
        output_fp, 
        index=None, 
        compression="gzip"
    )

print("Done")
"""
}


// Parse each individual alignment file and publish the final results
process collectFinalResults {
    
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: "copy", overwrite: true

    input:
        file "input/*"
    
    output:
        file "${params.output_prefix}.csv.gz"
    
"""
#!/usr/bin/env python3
import gzip
import json
import os
import pandas as pd
from shutil import copyfile

# Parse the list of files to join
csv_gz_list = [
    os.path.join("input", fp)
    for fp in os.listdir("input")
]
assert len(csv_gz_list) > 0

# Make sure all inputs are .csv.gz
for fp in csv_gz_list:
    assert fp.endswith(".csv.gz"), fp

# The output table will be named for one of the inputs
output_fp = "${params.output_prefix}.csv.gz"

# If there is only one file, just copy it to the output
if len(csv_gz_list) == 1:
    copyfile(
        csv_gz_list[0],
        output_fp
    )
else:

    print("Making a single output table")
    df = pd.concat([
        pd.read_csv(fp)
        for fp in csv_gz_list
    ], sort=True)

    print("Writing out to %s" % output_fp)
    df.to_csv(output_fp, index=None, compression="gzip")

print("Done")
"""
}

// Validate the input FASTA file format
process validateFASTA {
    container "${params.container__biopython}"
    label 'io_limited'

    input:
        tuple val(uuid), val(genome_name), path(fasta_gz)
    
    output:
        tuple val(uuid), val(genome_name), path("${uuid}.validated.fasta.gz")
    
"""
#!/usr/bin/env python3

import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Define the file paths, input and output
fp_in = "${fasta_gz}"
fp_out = "${uuid}.validated.fasta.gz"

# Define a counter for each unique header sequence
counter = {}

# Open the input and output files
with gzip.open(fp_in, "rt") as i, gzip.open(fp_out, "wt") as o:

    # Iterate over each input sequence
    for header, seq in SimpleFastaParser(i):

        # Truncate the header
        if len(header) > 20:
            header = header[:20]

        # Get the counter value for this header string
        i = counter.get(header, 0)

        # If the counter is greater than 0, add it to the output
        if i > 0:
            o.write(">%s.%d\\n%s\\n" % (header, i, seq))
        # Otherwise
        else:
            # Just write the header and sequence
            o.write(">%s\\n%s\\n" % (header, seq))

        print("Wrote out %d nucleotides for %s" % (len(seq), header))

        # Increment the counter
        counter[header] = i + 1

print("Done")

"""
}

// Make a results summary PDF
process summaryPDF {
    tag "Process final results"
    container "${params.container__plotting}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: "copy", overwrite: true

    input:
        file results_csv_gz
    
    output:
        file "${params.output_prefix}.pdf"
    
"""
#!/bin/bash

set -e

make_summary_figures.py "${results_csv_gz}" "${params.output_prefix}.pdf"
"""
}

// Annotate a genome with Prokka
process prokka {
    container "staphb/prokka:1.14.0-cv2"
    label "mem_medium"

    input:
    tuple val(genome_id), val(genome_name), file(fasta)

    output:
    tuple val(genome_id), val(genome_name), file("OUTPUT/${genome_id}.gbk.gz")

"""#!/bin/bash

set -euxo pipefail

echo Decompressing input file
gunzip -c "${fasta}" > INPUT.fasta

echo Running Prokka

prokka \
    --outdir OUTPUT \
    --prefix "${genome_id}" \
    --cpus ${task.cpus} \
    INPUT.fasta

echo Compressing outputs

gzip OUTPUT/*

echo Done
"""

}


// Extract the regions of each GBK which contains a hit
process extractGBK {
    container "${params.container__biopython}"
    label 'io_limited'
    publishDir params.output_folder, mode: 'copy', overwrite: true

    input:
        tuple val(genome_id), val(operon_context), val(operon_ix), val(contig_name), val(genome_name), file(annotation_gbk)
        each file(summary_csv)
    
    output:
        tuple val(operon_context), path("*/gbk/*gbk")
        path "*/faa/*faa.gz", optional: true

    script:
        template 'extractGBK.py'

}

// Make an interactive visual display for each operon context
process clinker {
    container "${params.container__clinker}"
    label 'mem_medium'
    publishDir "${params.output_folder}/html/", mode: 'copy', overwrite: true

    input:
        tuple val(operon_context), file(input_gbk_files)
    
    output:
        path "*html"

"""#!/bin/bash

set -Eeuxo pipefail

OUTPUT=\$(echo "${operon_context.replaceAll(/ :: /, '_')}" | sed 's/ (+)/_FWD/g' | sed 's/ (-)/_REV/g' | cut -c 1-200)

echo \$OUTPUT

ls -lahtr

clinker -p \$OUTPUT.html *gbk

"""


}


// Fetch genomes via FTP
process fetchFTP {
    tag "Download NCBI genomes by FTP"
    container 'quay.io/fhcrc-microbiome/wget:latest'
    label 'io_limited'
    maxForks params.ftp_threads

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

process linclust {
    tag "Cluster genes with similar sequences"
    container "${params.container__mmseqs}"
    label 'mem_medium'
    publishDir "${params.output_folder}/clusters/", mode: 'copy', overwrite: true
    
    input:
    file "input.genes.*.fasta.gz"
    
    output:
    file "clustered.genes.fasta.gz"
    file "clustered.genes.assignments.tsv.gz"
    file "clustered.genes.counts.tsv.gz"
    
"""
#!/bin/bash

set -Eeuo pipefail

# Combine input files
echo "Combining input files"
cat input.genes.* > input.genes.fasta.gz

# Make the MMSeqs2 database
echo "Running linclust"
mmseqs createdb input.genes.fasta.gz db

# Cluster the protein sequences
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.cluster_identity / 100} \
    -c ${params.cluster_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes clustered.genes.fasta --use-fasta-header

# Get the assignment of genes to clusters
mmseqs createtsv db db cluster_db clustered.genes.assignments.tsv

# Count up the number of times each representative was found
cat clustered.genes.assignments.tsv | cut -f 1 | sort | uniq -c | sort -nrk1 > clustered.genes.counts.tsv

echo "Compressing"
gzip clustered.genes.fasta
gzip clustered.genes.assignments.tsv
gzip clustered.genes.counts.tsv

echo "Done"

"""
}

process mashtree {
    container "${params.container__mashtree}"
    label 'mem_medium'
    publishDir "${params.output_folder}/tree/", mode: 'copy', overwrite: true
    
    input:
    file "*"
    
    output:
    file "mashtree.dnd"
    file "mashtree.tsv"
    
"""
#!/bin/bash

set -Eeuo pipefail

# Decompress the genome FASTA files
for f in *validated.fasta.gz; do
    gunzip -c \$f > \${f%.gz}
done

# Run mashtree
mashtree --outmatrix mashtree.tsv --numcpus ${task.cpus} *validated.fasta > mashtree.dnd

"""
}