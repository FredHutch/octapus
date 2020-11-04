#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
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

# Determine whether the operon is on the reverse strand
# Start by parsing the operon context to see what strand(s) 
# are expected for each gene
target_gene_strand = defaultdict(list)
for field in '${operon_context}'.split(" :: "):
    gene_name, strand = field.split(" ", 1)
    assert strand in ["(-)", "(+)"]
    target_gene_strand[gene_name].append(
        strand[1] # Remove the parentheses
    )
# Now for the genes in the operon on this contig, see
# if they are consistent either with + or - overall orientation
fwd_score = 0
rev_score = 0
for _, r in operon_df.iterrows():
    assert r["strand"] in ["+", "-"], r
    if r["strand"] in target_gene_strand[r["gene_name"]]:
        fwd_score += 1
    if {"+": "-", "-": "+"}.get(r["strand"]) in target_gene_strand[r["gene_name"]]:
        rev_score += 1
assert fwd_score > 0 or rev_score > 0, operon_df.reindex(columns=["gene_name", "strand"])
if fwd_score >= rev_score:
    print("Operon is on the forward strand")
else:
    print("Operon is on the reverse strand, so the output GBK will be reverse complemented")
    recs = [
        r.reverse_complement()
        for r in recs
    ]

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
