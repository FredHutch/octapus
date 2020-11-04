#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd

print("Processing {genome_id} -- {genome_name}")

print("Reading in annotations from {annotation_gbk}")
recs = [rec for rec in SeqIO.parse('{annotation_gbk}', "genbank")]
print(f"Read in {len(recs):,} records")

print("Reading in operon data from {summary_csv}")
operon_df = pd.read_csv("{summary_csv}")
print(f"Read in operon information for {operon_df.shape[0]:,} genes")
