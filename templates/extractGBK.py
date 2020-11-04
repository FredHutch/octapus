#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd

# Capture input variables populated by the Nextflow controller
genome_id = '${genome_id}'
operon_context = '${operon_context}'
contig_name = '${contig_name}'
genome_name = '${genome_name}'
annotation_gbk = '${annotation_gbk}'
summary_csv = '${summary_csv}'

print(f"Processing {genome_id} -- {genome_name}")

print(f"Reading in annotations from {annotation_gbk}")
recs = [rec for rec in SeqIO.parse(annotation_gbk, "genbank")]
print(f"Read in {len(recs):,} records")

print(f"Reading in operon data from {summary_csv}")
operon_df = pd.read_csv(annotation_gbk)
print(f"Read in operon information for {operon_df.shape[0]:,} genes")
