#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import gzip

print("Processing ${genome_id} -- ${genome_name}")

print("Reading in annotations from ${annotation_gbk}")
recs = [
    rec
    for rec in SeqIO.parse(
        gzip.open(
            '${annotation_gbk}', 
            'rt'
        ), 
        "genbank"
    )
]
print("Read in %d records" % len(recs))

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
)

print(operon_df)
