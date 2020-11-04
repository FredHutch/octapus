#!/usr/bin/env python3

from Bio import SeqIO
import gzip
import os
import pandas as pd

print("Processing ${genome_id} -- ${genome_name}")
print("Contig: ${contig_name}")
print("Operon Context: ${operon_context}")
print("Operon index: ${operon_ix}")

print("Reading in operon data from ${summary_csv}")
operon_df = pd.read_csv("${summary_csv}")
print("Read in operon information for %d genes" % operon_df.shape[0])

print("Filtering down to operons from this genome and contig")
print(operon_df.head())
operon_df = operon_df.query(
    "genome_id == '${genome_id}'"
).query(
    "contig_name == '${contig_name}'"
).query(
    "operon_context == '${operon_context}'"
).query(
    "operon_ix == '${operon_ix}'"
)

# Get the smallest and largest coordinate
start_pos = operon_df.reindex(columns=["contig_start", "contig_end"]).min(axis=1).min()
end_pos = operon_df.reindex(columns=["contig_start", "contig_end"]).max(axis=1).max()
print("Operon spans %d to %d on ${contig_name}" % (start_pos, end_pos))

# Add the window on each side
start_pos = start_pos - ${params.annotation_window}
if start_pos < 1:
    start_pos = 1

end_pos = end_pos + ${params.annotation_window}

print("Start position: %d" % start_pos)
print("End position: %d" % end_pos)

print("Filtering annotations from ${annotation_gbk}")
recs = [
    rec[start_pos:min(end_pos, len(rec.seq) - 1)]
    for rec in SeqIO.parse(
        gzip.open(
            '${annotation_gbk}',
            'rt'
        ),
        "genbank"
    )
    if rec.id == '${contig_name}'
]
print("Read in %d records" % len(recs))

# Format the folder name
operon_folder = "${operon_context}".replace(
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
os.mkdir(operon_folder)
os.mkdir(os.path.join(operon_folder, "gbk"))

# Define the output file name
output_fp = os.path.join(
    operon_folder,
    "gbk",
    "${genome_id}-${contig_name}-${operon_ix}.gbk"
)
print("Writing out to %s" % output_fp)
with open(output_fp, "wt") as handle:
    SeqIO.write(
        recs,
        handle,
        "genbank"
    )
