// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"


// Parse each individual alignment file
process collectResults {
    
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

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
    
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    errorStrategy "retry"

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