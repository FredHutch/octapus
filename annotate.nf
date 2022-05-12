#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.octapus_results = false
params.annotation_window = 10000
params.ftp_threads = 100

// Import modules
include {
    validateFASTA;
    prokka;
    extractGBK;
    clinker;
    sanitize_manifest;
    fetchFTP;
} from './modules/modules' addParams(
    output_folder: params.output_folder,
    ftp_threads: params.ftp_threads,
    annotation_window: params.annotation_window,
)

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/OCTAPUS/annotate.nf <ARGUMENTS>
    
    Required Arguments:
      --genomes             CSV file listing genomes (from https://www.ncbi.nlm.nih.gov/genome/browse)
      --octapus_results       Spreadsheet containing OCTAPUS outputs (in CSV format)
      --output_folder       Folder to write output files to

    Optional Arguments:
      --annotation_window   The additional area on either side of the operon to annotate (in bp) (default: 10000)

    
    Given a set of OCTAPUS results (or a subset of just those results the user is particularly interested in),
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
    if (params.help || !params.genomes || !params.output_folder || !params.octapus_results){
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
    if (file(params.octapus_results).isEmpty()){
        log.info"""
        The specified --octapus_results file cannot be found!
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

    // Process all FASTA inputs to make sure that their format is valid
    validateFASTA(
        joined_fasta_ch
    )

    // Parse the OCTAPUS outputs provided by the user
    parse_spreadsheet(
        Channel.fromPath(
            params.octapus_results
        )
    )

    // Make a channel with the genomes containing all of the operons identified in this analysis
    genome_ch = parse_spreadsheet.out.map {
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
    annotation_ch = parse_spreadsheet.out.map {
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
        parse_spreadsheet.out
    )

    // Make a clinker webpage for each operon context
    clinker(
        extractGBK.out[0].groupTuple().filter({
            it -> it[1].size() > 1
        })
    )

}

process parse_spreadsheet {
    container "${params.container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file input_spreadsheet
    
    output:
        file "${input_spreadsheet}.csv"

"""#!/usr/bin/env python3

import pandas as pd

# Parse the filepath
input_spreadsheet = '${input_spreadsheet}'

# Use the file ending to determine how it is read

# If the input is a gzip-compressed CSV
if input_spreadsheet.endswith(".csv.gz"):
    df = pd.read_csv(
        input_spreadsheet,
        sep=",",
        compression="gzip"
    )

# If the input is a uncompressed CSV
elif input_spreadsheet.endswith(".csv"):
    df = pd.read_csv(
        input_spreadsheet,
        sep=",",
    )

# If the input is an Excel spreadsheet
elif input_spreadsheet.endswith(".xlsx"):

    df = pd.read_excel(
        input_spreadsheet,
    )

else:

    assert False, "Could not parse file ending: %s" % input_spreadsheet

# Write out to a file
df.to_csv(
    "${input_spreadsheet}.csv",
    index=None
)

"""

}