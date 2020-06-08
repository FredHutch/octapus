#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.output_prefix = false
params.csv_list = false

// Import modules
include summaryPDF from './modules/modules' params(
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

    nextflow run FredHutch/BOFFO/join_results.nf <ARGUMENTS>
    
    Required Arguments:
      --csv_list            Comma-separated list of CSV results files to join
      --output_folder       Folder to write output files to
      --output_prefix       Prefix to use for output file names

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || !params.output_folder || !params.output_prefix || !params.csv_list){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // Join together the CSV files
    joinCSVs(
        Channel.of(
            params.csv_list.split(",")
        ).map {
            r -> file(r)
        }.toSortedList()
    )

    // Make the summary PDF
    summaryPDF(
        joinCSVs.out
    )

}

// #############
// # PROCESSES #
// #############


// Join the list of CSV input files
process joinCSVs {
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

    input:
        path csv_list
    
    output:
        path "${params.output_prefix}.joined.csv.gz"
    
"""#!/usr/bin/env python3

import pandas as pd
import os

csv_list = "${csv_list}".split(" ")

for fp in csv_list:
    print("Checking for input: %s" % fp)
    assert os.path.exists(fp)

df = pd.concat([
    pd.read_csv(fp)
    for fp in csv_list
])


df.to_csv("${params.output_prefix}.joined.csv.gz", index=None, sep=",")
"""
}
