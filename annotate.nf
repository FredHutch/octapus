#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.boffo_results = false
params.genomes = false
params.annotation_window = 10000

// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__biopython = "quay.io/fhcrc-microbiome/biopython-pandas:latest"
container__plotting = "quay.io/fhcrc-microbiome/boffo-plotting:latest"
container__clinker = "quay.io/fhcrc-microbiome/clinker:v0.0.4"

// Import modules
include {
    prokka;
    extractGBK;
    clinker;
} from './modules/modules' params(
    output_prefix: params.output_prefix,
    output_folder: params.output_folder,
    annotation_window: params.annotation_window,
    container__biopython: container__biopython,
    container__clinker: container__clinker,
)

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/BOFFO/annotate.nf <ARGUMENTS>
    
    Required Arguments:
      --genomes             CSV file listing genomes (from https://www.ncbi.nlm.nih.gov/genome/browse)
      --boffo_results       Spreadsheet containing BOFFO outputs (in CSV format)
      --output_folder       Folder to write output files to

    Optional Arguments:
      --annotation_window   The additional area on either side of the operon to annotate (in bp) (default: 10000)

    
    Given a set of BOFFO results (or a subset of just those results the user is particularly interested in),
    this utility will extract the annotations from the surrounding region and generate an interactive
    visualization displaying the position and orientation of those genes, including ones which were
    not included in the original search.


    Citations:

    Seemann T. Prokka: rapid prokaryotic genome annotation
    Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063

    Gilchrist, C.L.M, 2020. clinker: Easy gene cluster comparison figure generator.

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || !params.genomes || !params.output_folder || !params.boffo_results){
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

    // Make sure that the input files exist
    if (file(params.boffo_results).isEmpty()){
        log.info"""
        The specified --boffo_results file cannot be found!
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

    // Make a new channel joining the FTP files with local files
    joined_fasta_ch = fetchFTP.out.mix(
        split_genome_ch.local.map {
            r -> [
                r["uri"].split("/")[-1],
                r["#Organism Name"],
                file(r["uri"])
            ]
        }
    )

    // Parse the BOFFO outputs provided by the user
    boffo_results_ch = Channel.fromPath(
        params.boffo_results
    )

    // Make a channel with the genomes containing all of the operons identified in this analysis
    genome_ch = boffo_results_ch.map {
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
        joined_fasta_ch
    )

    // Annotate these genomes with prokka
    prokka(
        genome_ch
    )

    // Make a channel with each unique operon per genome and contig
    annotation_ch = boffo_results_ch.map {
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
    ).join( // Add the genomes
        prokka.out
    )

    // Extract the regions of the GBK files which operons fall into
    extractGBK(
        annotation_ch,
        boffo_results_ch
    )

    // Make a clinker webpage for each operon context
    clinker(
        extractGBK.out.groupTuple()
    )

}