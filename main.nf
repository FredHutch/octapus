#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.output_prefix = false
params.genomes = false
params.genomes_local = false
params.operon = false
params.operon_list = false
params.min_identity = 90
params.min_coverage = 50
params.max_operon_gap = 10000
params.batchsize = 100
params.ftp_threads = 100
params.max_evalue = 0.001
params.annotations = false
params.annotation_window = 10000
params.cluster_identity = 80
params.cluster_coverage = 50
params.mashtree = false

// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__biopython = "quay.io/fhcrc-microbiome/biopython-pandas:latest"
container__plotting = "quay.io/fhcrc-microbiome/boffo-plotting:latest"
container__clinker = "quay.io/fhcrc-microbiome/clinker:v0.0.16--1"
container__mmseqs = "quay.io/fhcrc-microbiome/mmseqs2:version-12"
container__mashtree = "quay.io/hdc-workflows/mashtree:1.2.0"

// Import modules
include {
    collectResults as collectResultsRound1;
    collectResults as collectResultsRound2;
    collectFinalResults;
    validateFASTA;
    summaryPDF;
    prokka;
    extractGBK;
    linclust;
    clinker;
    sanitize_manifest;
    sanitize_manifest_local;
    mashtree;
} from './modules/modules' params(
    output_prefix: params.output_prefix,
    output_folder: params.output_folder,
    ftp_threads: params.ftp_threads,
    annotation_window: params.annotation_window,
    cluster_identity: params.cluster_identity,
    cluster_coverage: params.cluster_coverage,
    container__pandas: container__pandas,
    container__plotting: container__plotting,
    container__biopython: container__biopython,
    container__clinker: container__clinker,
    container__mmseqs: container__mmseqs,
    container__mashtree: container__mashtree,
)

// Import modules fetchFTP for genome FASTA files
include {
    fetchFTP;
} from './modules/modules' params(
    ftp_threads: params.ftp_threads,
    ftp_suffix: "_genomic.fna.gz",
    local_suffix: ".fasta.gz",
)

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/OCTAPUS <ARGUMENTS>
    
    Required Arguments:
    At least ONE of:
      --genomes             CSV file listing genomes (from https://www.ncbi.nlm.nih.gov/genome/browse)
      --genomes_local       CSV file listing local genomes, with columns '#Organism name' and 'uri'
    And:
      --operon              Amino acid sequences to search for, in multi-FASTP format
        OR
      --operon_list         Flag to input multiple sequences per gene -- see description below    
    And:      
      --output_folder       Folder to write output files to
      --output_prefix       Prefix to use for output file names

    Optional Arguments:
      --min_identity        Percent identity threshold used for alignment (default: 90)
      --min_coverage        Percent coverage threshold used for alignment (default: 50)
      --max_operon_gap      Maximum gap between genes in the same 'operon' (only used for the 'operon_context' output column) (default: 10000)
      --batchsize           Number of samples to join in each batch (default: 100)
      --max_evalue          Maximum E-value threshold used to filter initial alignments (default: 0.001)
      --mashtree            If specified, generate a tree (via mashtree) of all genomes (default: do not run)
      --annotations         If specified, annotate the regions of all genomes which contain operons
      --annotation_window   The additional area on either side of the operon to annotate (in bp) (default: 10000)
      --cluster_identity    Percent similarity used to cluster adjacent genes (default: 80)
      --cluster_coverage    Percent coverage used to cluster adjacent genes (default: 50)

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

    As an alternative approach, we suggest using the entrypoint 'FredHutch/OCTAPUS/annotate.nf',
    which allows the user to annotate just a subset of the results from a set of OCTAPUS alignments.


    Citations:

    Seemann T. Prokka: rapid prokaryotic genome annotation
    Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063

    Gilchrist, C.L.M, 2020. clinker: Easy gene cluster comparison figure generator.

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help  || !params.output_folder || !params.output_prefix){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    // The user must specify genomes and/or genomes_local
    if (
        (!params.genomes && !params.genomes_local) 
    ){
        log.info"""

        Please specify --genomes and/or --genomes_local that are not empty.

        Use the --help flag for more usage information.
        
        """.stripIndent()

        // Exit out and do not run anything else
        exit 0        

    }
    
    // The user must specify _either_ --genomes or --operon_list
    if (!params.operon_list && !params.operon){
        log.info"""

        Please specify either --operon or --operon_list.

        Use the --help flag for more usage information.
        
        """.stripIndent()

        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify _either_ --operon or --operon_list, but not both
    if (params.operon_list && params.operon){
        log.info"""

        Please specify either --operon or --operon_list, but not both.

        Use the --help flag for more usage information.

        """.stripIndent()

        // Exit out and do not run anything else
        exit 0
    }

    if ( "${params.genomes}" != "false" ) {
        // Parse the manifest to get a name and FTP prefix for each genome
        sanitize_manifest(
            file("${params.genomes}")
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

        // Add both sets of files to the combined channel
        fetchFTP
            .out
            .mix(
                split_genome_ch.local.map {
                    r -> [
                        r["uri"].split("/")[-1],
                        r["#Organism Name"],
                        file(r["uri"])
                    ]
                }
            )
            .set{genomes_ch}

    } else {

        genomes_ch = Channel.empty()

    }

    if ( "${params.genomes_local}" != "false" ) {
        // Parse the manifest to get a name and FTP prefix for each genome
        sanitize_manifest_local(
            file("${params.genomes}")
        )
        sanitize_manifest_local.out.map {
            r -> r.splitCsv(
                header: true
            )
        }.flatten(
        ).map {
            r -> [
                    r["uri"].split("/")[-1],
                    r["#Organism Name"],
                    file(r["uri"])
                ]
        }
        .set(genomes_local_ch)

    } else {

        genomes_local_ch = Channel.empty()

    }

    genomes_ch
        .mix(genomes_local_ch)
        .set{joined_fasta_ch}

    // Process all FASTA inputs to make sure that their format is valid
    validateFASTA(
        joined_fasta_ch
            .ifEmpty { error "No FASTA inputs found" }
    )
    if (params.mashtree) {
        // Create a tree summarizing whole-genome similarity
        mashtree(
            validateFASTA.out.map({
                it -> it[2]
            }).toSortedList()
        )
    }

    // Each gene in the operon is represented by a single sequence in a
    // multi-FASTA which contains all of the genes in the operon
    if (params.operon){

        // Point to the operon file
        operon_fasta = file("${params.genomes}")    // Align the operon against each genome

        // Make sure the file is not empty
        if (operon_fasta.isEmpty()){
            log.info"""
            The specified --operon file cannot be found!
            """.stripIndent()
            exit 1
        }

        // Simply run BLAST on each of the genomes
        runBLAST(
            validateFASTA.out,
            operon_fasta
        )

        // Parse each individual alignment
        parseAlignments(
            runBLAST.out
        )

    }else{

        // Instead, each gene is specified as its own multi-FASTA
        makePSSM(
            Channel.from(
                params.operon_list.split(",")
            ).flatten(
            ).map {
                r -> [
                    r.split("=")[0], file(r.split("=")[1])
                ]
            }
        )

        // Run PSIBLAST for each gene against this genome
        runPSIBLAST(
            validateFASTA.out,
            makePSSM.out[0].toSortedList()
        )

        // Parse each individual alignment
        parseAlignments(
            runPSIBLAST.out
        )

    }
    
    // Join the parsed alignments with the FASTA for each genome
    // For each alignment, extract the sequence of the aligned region
    extractAlignments(
        validateFASTA.out.join(
            parseAlignments.out.map {
                r -> [r.name.replaceAll(/.csv.gz/, ""), r]
            }
        )
    )

    // Collect results in rounds
    collectResultsRound1(
        extractAlignments
            .out
            .ifEmpty{ error "No alignments found" }
            .collate(params.batchsize)
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

    // If the --annotations flag is set
    if (params.annotations){

        // Make a channel with the genomes containing all of the operons identified in this analysis
        genome_ch = collectFinalResults.out.map {
            r -> r.splitCsv(
                header: true
            )
        }.flatten(
        ).map { // Just keep the genome ID across all hits
            r -> [
                r["genome_id"],
            ]
        }.unique( // Drop duplicate entries
        ).join( // Add the genomes
            validateFASTA.out
        )

        // Annotate these genomes with prokka
        prokka(
            genome_ch
        )

        // Make a channel with each unique operon per genome and contig
        // Combine the operon coordinates with the genome files
        annotation_ch = prokka.out.cross(
            collectFinalResults.out.map {
                r -> r.splitCsv(
                    header: true
                )
            }.flatten(
            ).map { // Just keep the genome ID across all hits
                r -> [
                    r["genome_id"],
                    r["operon_context"],
                    r["operon_ix"],
                    r["contig_name"],
                ]
            }.unique( // Drop duplicate entries
            )
        ).map {
            i -> [
                i[0][0],
                i[1][1],
                i[1][2],
                i[1][3],
                i[0][1],
                i[0][2]
            ]
        }

        // Extract the regions of the GBK files which operons fall into
        extractGBK(
            annotation_ch,
            collectFinalResults.out
        )

        // Make a clinker webpage for each operon context,
        // filtering to those operon contexts which contain > 1 representative
        clinker(
            extractGBK.out[0].groupTuple(
            ).filter({
                it -> it[1].size() > 1
            })
        )

        // Cluster all of the adjacent genes by amino acid similarity
        linclust(
            extractGBK.out[1].flatten().toSortedList()
        )
        
    }

}

// #############
// # PROCESSES #
// #############


// Align the operon against each individual genome
process runBLAST {
    tag "Align operon"
    container 'quay.io/fhcrc-microbiome/blast:latest'
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
    -evalue ${params.max_evalue}

echo "Compressing alignment file"

gzip ${uuid}.aln

echo Done
"""
}

// Make a PSSM for each gene in the list of queries
process makePSSM {
    tag "Identify conserved positions"
    container 'quay.io/fhcrc-microbiome/blast:latest'
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true

    input:
        tuple val(gene_name), path(fasta)
    
    output:
        file "pssm/${gene_name}.pssm"
        file "pssm/${gene_name}.internal.aln"
    
"""
#!/bin/bash
set -e

echo "Building PSSM for ${gene_name}"

# First, extract the first gene in the list
cat ${fasta} \
| tr '\n' '\t' \
| tr '>' '\n' \
| head -2 \
| tail -n 1 \
| sed 's/^/>/' \
| tr '\t' '\n' \
| sed 's/>.*/>${gene_name}/' \
> first_sequence.fasta

# Second, extract everything after the first gene in the list
cat ${fasta} \
| tr '\n' '\t' \
| tr '>' '\n' \
| awk 'NR > 2' \
| sed 's/^/>/' \
| tr '\t' '\n' \
> other_sequences.fasta

# Third, run PSIBLAST to generate the PSSM
mkdir pssm
psiblast \
    -query first_sequence.fasta \
    -subject other_sequences.fasta \
    -out_pssm pssm/${gene_name}.pssm \
    -save_pssm_after_last_round \
    -out pssm/${gene_name}.internal.aln

echo Done
"""
}

// Align the PSSM for each gene in the operon against each individual genome
process runPSIBLAST {
    tag "Align PSSM for operon"
    container 'quay.io/fhcrc-microbiome/blast:latest'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(uuid), val(genome_name), path(fasta_gz)
        file pssm_list
    
    output:
        tuple val(uuid), val(genome_name), file("${uuid}.aln.gz")
    
"""
#!/bin/bash
set -e

echo "Running alignment for ${uuid}"

ls -lahtr

for PSSM in *pssm; do

    echo "Processing \$PSSM"

    tblastn \
        -in_pssm \$PSSM \
        -subject <(gunzip -c ${fasta_gz}) \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen" \
        -evalue ${params.max_evalue} \
        >> ${uuid}.aln

done

echo "Compressing alignment file"

gzip ${uuid}.aln

echo Done
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

    # Use the lexographically lower string
    return min(fwd_str, rev_str)

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
    for operon_ix, operon in enumerate(operons):
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

        # Assign the lexographically lower value for the operon structure string
        for _ in operon:
            output.append(dict(
                operon_ix=operon_ix,
                operon_context=min(fwd_str, rev_str),
            ))
            
    return pd.DataFrame(output, index=df.index)

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
    df = pd.concat([
        df,
        operon_context(df, ${params.max_operon_gap})
    ], axis=1).assign(
        genome_context = genome_context(df)
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
            "operon_ix",
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
        path "*/seq/*.gz"
    
"""
#!/usr/bin/env python3

from collections import defaultdict
import gzip
import os
import pandas as pd

# We must limit the length of each file name to 200 characters
folder_prefix_max_len = 200

# Read in all of the results for this dataset
df = pd.read_csv("${results_csv_gz}")

# We must limit the length of each folder name to 200 characters
folder_prefix_counter = defaultdict(int)

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

    # If the file prefix is longer than the maximum length
    if len(operon_folder) > folder_prefix_max_len:

        # Truncate down to the maximum length
        operon_folder = operon_folder[:folder_prefix_max_len]

    # Add the file prefix to the counter
    folder_prefix_counter[operon_folder] += 1

    # If there is more than one of these file prefixes
    if folder_prefix_counter[operon_folder] > 1:

        # Then add the counter to the file prefix
        operon_folder = "%s_%s" % (operon_folder, folder_prefix_counter[operon_folder])

    print("Writing out %d sequences for %s" % (operon_df.shape[0], operon_structure))
    os.mkdir(operon_folder)
    os.mkdir(os.path.join(operon_folder, "seq"))

    # For each gene, make a FASTA file with the nucleotide sequence
    for gene_name, gene_df in df.groupby("gene_name"):

        # Format the file name prefix
        file_prefix = "%s/seq/%s.alignments" % (operon_folder, gene_name)

        # Write out the nucleotide sequences
        with gzip.open("%s.fasta.gz" % file_prefix, "wt") as fo:
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
        with gzip.open("%s.fastp.gz" % file_prefix, "wt") as fo:
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
