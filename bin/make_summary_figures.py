#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages

sns.set_style("whitegrid")

def calc_indiv_operon(operon_df):

    # Pivot the operon to be oriented by contig position
    contig_df = operon_df.pivot_table(
        index="contig_name",
        columns="gene_name",
        values=["contig_left", "contig_right", "contig_span", "strand"],
        aggfunc=max
    )

    # Aggregate some information about the gene_order
    gene_order = pd.DataFrame([
        {
            "gene_name": n.replace(" (+)", "").replace(" (-)", ""),
            "gene_ix": ix,
            "strand": "+" if " (+)" in n else "-"
        }
        for ix, n in enumerate(operon_df.name.split(" :: "))
    ]).set_index("gene_ix")

    # Set the upstream and downstream genes
    gene_order = gene_order.assign(
        upstream=pd.Series({
            ix: gene_order.loc[ix - 1, "gene_name"] if ix > 0 else None
            for ix in range(gene_order.shape[0])
        }),
        downstream=pd.Series({
            ix: gene_order.loc[ix + 1,
                               "gene_name"] if ix < gene_order.shape[0] - 1 else None
            for ix in range(gene_order.shape[0])
        })
    )

    # Reverse the coordinates for contigs with the opposite orientation
    contig_df = contig_df.apply(
        lambda r: (r * -1).rename(index={"contig_left": "contig_right", "contig_right": "contig_left"}
                                  ) if r["strand"][gene_order.gene_name.values[0]] != gene_order.strand.values[0] else r,
        axis=1
    ).drop(columns="strand")

    # Make sure that all of the genes are present on the contig
    for gene_name in gene_order["gene_name"].values:
        if gene_name not in contig_df["contig_left"].columns.values:
            return
        if gene_name not in contig_df["contig_right"].columns.values:
            return

    df = pd.DataFrame([
        {
            "gene_name": r["gene_name"],
            "gene_ix": gene_ix,
            "gene_length": contig_df["contig_span"][r["gene_name"]].abs().median(),
            "padding_left": None if gene_ix == 0 else (contig_df["contig_left"][r["gene_name"]] - contig_df["contig_right"][r["upstream"]]).median(),
            "padding_right": None if gene_ix == gene_order.shape[0] - 1 else (contig_df["contig_left"][r["downstream"]] - contig_df["contig_right"][r["gene_name"]]).median(),
            "n_genomes": operon_df["genome_id"].unique().shape[0]
        }
        for gene_ix, r in gene_order.iterrows()
    ])

    return df


def average_operon_structure(csv_fp):
    df = pd.read_csv(csv_fp)

    # Calculate the 'left' and 'right' for each gene
    df = df.assign(
        contig_left=df[["contig_start", "contig_end"]].min(axis=1),
        contig_right=df[["contig_start", "contig_end"]].max(axis=1),
        contig_span=(df["contig_start"] - df["contig_end"]).abs()
    )

    return df.groupby(
        "operon_context"
    ).apply(
        calc_indiv_operon
    ).reset_index(
    ).drop(
        columns="level_1"
    )

def draw_arrow(
    x_pos,
    body_width,
    head_width=0.3,
    head_height=2,
    body_height=1,
    **kwargs
):

    # Start from the bottom left and go clockwise
    points = [
        [x_pos, -body_height / 2],
        [x_pos, body_height / 2],
        [x_pos + (body_width * (1 - head_width)), body_height / 2],
        [x_pos + (body_width * (1 - head_width)), head_height / 2],
        [x_pos + body_width, 0],
        [x_pos + (body_width * (1 - head_width)), -head_height / 2],
        [x_pos + (body_width * (1 - head_width)), -body_height / 2],
        [x_pos, -body_height / 2]
    ]

    return patches.Polygon(
        points,
        **kwargs
    )


def plot_single_operon(
    operon_name,
    operon_df,
    cmap,
    x_padding=0.1,
    y_prop=0.2,
    pdf=None
):
    # Sort the table
    operon_df.sort_values(
        by="gene_ix",
        inplace=True
    )
    # Set the index
    operon_df.set_index(
        "gene_ix",
        inplace=True
    )

    # Set up the figure
    fig, ax = plt.subplots()

    # Keep track of the horizontal position
    cumulative_x_position = 0

    # Iterate over each gene
    for ix, r in operon_df.iterrows():
        # Move the X position
        if ix > 0:
            cumulative_x_position += r["padding_left"]

        # Draw an arrow
        rect = draw_arrow(
            cumulative_x_position,
            r["gene_length"],
            facecolor=cmap[r["gene_name"]],
            label=r["gene_name"],
            head_width=0.2,
            head_height=2,
            body_height=1
        )

        # Add the patch to the Axes
        ax.add_patch(rect)

        # Move the X position
        cumulative_x_position += r["gene_length"]
        if ix < operon_df.shape[0] - 1:
            cumulative_x_position += r["padding_right"]

    plt.xlim(-x_padding * cumulative_x_position,
             (1 + x_padding) * cumulative_x_position)
    plt.ylim(-(0.5 / y_prop), (0.5 / y_prop))
    ax.yaxis.set_visible(False)
    plt.xlabel("Operon position (bps)")
    plt.legend()
    plt.title(operon_name)
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    plt.show()


def plot_all_operons(
    all_operon_df,
    pdf=None
):
    genome_counts = all_operon_df.reindex(columns=[
        "operon_context",
        "n_genomes"
    ]).drop_duplicates(
    ).set_index(
        "operon_context"
    )[
        "n_genomes"
    ].sort_values(
        ascending=True
    )

    genome_counts.plot(kind="barh")
    plt.xlabel("Number of genomes")
    plt.ylabel("")
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    plt.show()

    genome_counts.apply(np.log10).plot(kind="barh")
    plt.xlabel("Number of genomes (log10)")
    plt.ylabel("")
    if pdf is not None:
        pdf.savefig(bbox_inches="tight")
    plt.show()

    # Assign a color for each gene
    cmap = dict(zip(
        all_operon_df["gene_name"].unique(),
        sns.color_palette(
            "colorblind", 
            all_operon_df["gene_name"].unique().shape[0]
        )
    ))

    for operon_name in genome_counts.index.values[::-1]:
        plot_single_operon(
            operon_name,
            all_operon_df.query("operon_context == '{}'".format(operon_name)),
            cmap,
            pdf=pdf
        )


if __name__ == "__main__":
    # Get the input and output file paths
    input_fp = sys.argv[1]
    output_fp = sys.argv[2]

    assert os.path.exists(input_fp), "Cannot find input file %s" % input_fp
    assert os.path.exists(output_fp) is False, "Output file path %s already exists" % output_fp
    assert output_fp.endswith(".pdf"), "Output file path must end with '.pdf'"

    # Read in the operon structure
    operon_df = average_operon_structure(input_fp)

    with PdfPages(output_fp) as pdf:
        plot_all_operons(operon_df, pdf=pdf)
