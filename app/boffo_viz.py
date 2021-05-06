#!/usr/bin/env python3

import argparse
from Bio import Phylo
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State, MATCH, ALL
from dash.exceptions import PreventUpdate
from flask import Flask, send_file
from functools import lru_cache
import json
import logging
import os
import pandas as pd
from plotly.colors import sequential
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, leaves_list
from uuid import uuid4

##################
# SET UP LOGGING #
##################

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [BOFFO] %(message)s'
)
logger = logging.getLogger('BOFFO')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description='Visualize the results of the BOFFO analysis pipeline'
)

# Add the arguments
parser.add_argument(
    '--csv',
    type=str,
    required=True,
    help='BOFFO output (*.csv.gz)'
)
parser.add_argument(
    '--dnd',
    type=str,
    required=True,
    help='Whole-genome comparison tree (tree/*.dnd)'
)
parser.add_argument(
    '--tsv',
    type=str,
    required=True,
    help='Whole-genome distance matrix (tree/*.tsv)'
)
parser.add_argument(
    '--names',
    type=str,
    default=None,
    help='(Optional) File containing a name for each genome (CSV format, columns are "#Organism Name" and "uri")'
)

# Parse the arguments
args = parser.parse_args()


############################
# PARSE THE TREE STRUCTURE #
############################
class TreeLayout:
    """Object to parse tree layout from a Bio.Phylo object."""

    def __init__(self, tree):
        """Parse a Bio.Phylo object."""

        # Get the total number of termini in the tree
        n_termini = len(tree.clade.get_terminals())

        # Set up a dict with the X-Y and parent for each node
        self.nodes = dict(
            root=dict(
                x=0,
                y=n_termini / 2.,
                terminal=False,
                parent=None
            )
        )

        # Number the internal nodes starting with 0
        self.internal_n = 0

        # Iteratively parse the clades beneath the root
        self.parse_clade(
            tree.clade,
            parent='root',
            max_y=n_termini,
        )

    def parse_clade(self, clade, parent=None, x=0, max_y=1, min_y=0):
        """Add positions for the children of a clade."""

        # Get the total number of termini below this clade
        n_termini = len(clade.get_terminals())

        # Compute the amount of vertical space to allot each terminal
        terminal_span = (max_y - min_y) / n_termini

        # For each child below this clade
        for child_ix, child in enumerate(clade.clades):

            # If the node does not have a name
            if child.name is None:

                # Name it using the counter for internal nodes
                child_name = f"internal-{self.internal_n}"

                # Increment the counter
                self.internal_n += 1

            # If the node does have a name
            else:

                # Set the variable `child_name`
                child_name = child.name.replace(".validated", "")

            # Set the X position by adding the branch length
            child_x = x + child.branch_length

            # Set the Y position based on the number of termini
            child_span = len(child.get_terminals()) * terminal_span
            child_y = max_y - (child_span / 2.)

            # Add the node to the data
            self.nodes[
                child_name
            ] = dict(
                x=child_x,
                y=child_y,
                terminal=child.is_terminal(),
                parent=parent
            )

            # Add the children of this node
            self.parse_clade(
                child,
                parent=child_name,
                x=child_x,
                max_y=max_y,
                min_y=max_y - child_span
            )

            # Lower the ceiling of the remaining Y space available
            # to parse the next clade
            max_y = max_y - child_span


#################
# READ THE DATA #
#################

for fp in [args.csv, args.dnd, args.tsv]:
    assert os.path.exists(fp), f"File not found: {fp}"

# Read in the table
operons = pd.read_csv(args.csv)

# Get the human-readable name for each genome
genome_name_dict = operons.reindex(
    columns=["genome_id", "genome_name"]
).drop_duplicates(
)

# If there is a duplicate ID
if genome_name_dict["genome_id"].unique().shape[0] < genome_name_dict.shape[0]:

    # Then make a dummy empty dict
    genome_name_dict = dict()

# Otherwise
else:

    # Make a dict linking each genome ID to the human-readable name
    genome_name_dict = genome_name_dict.set_index(
        "genome_id"
    )[
        "genome_name"
    ].to_dict()

# If the user provided an additional names file
if args.names is not None:

    logger.info(f"User provided names in {args.names}")

    # And that file is present
    if os.path.exists(args.names):

        # Read in the table
        names_df = pd.read_csv(args.names)

        # Select the columns for #Organism Name and uri
        names_df = names_df.reindex(
            columns=["#Organism Name", "uri"]
        )

        # Iterate over each line
        for _, r in names_df.iterrows():

            # If there are missing entries for either column
            if pd.isnull(r["#Organism Name"]) or pd.isnull(r["uri"]):

                # Skip the row
                continue

            else:

                # Add this name to the dict
                genome_name_dict[
                    r["uri"].split("/")[-1]
                ] = r["#Organism Name"]

    else:
        logger.info(f"File not found: {args.names}")

# Get the list of all operons
operon_list = operons.assign(
    operon_len = operons["operon_context"].apply(
        lambda s: len(s.split(" "))
    )
).sort_values(
    by="operon_len",
    ascending=False
)[
    "operon_context"
].drop_duplicates(
).values

# Read in the tree
tree = Phylo.read(args.dnd, "newick")

# Perform midpoint rooting on the tree
tree.root_at_midpoint()

# Arrange tips with ladderize
tree.ladderize()

# Compute the tree layout
tree_layout = TreeLayout(tree)

# Read in the distance matrix
dist = pd.read_csv(args.tsv, sep="\t")

##################
# SET UP THE APP #
##################

FONT_AWESOME = {
    "src": "https://kit.fontawesome.com/f8b0dec9e6.js",
    "crossorigin": "anonymous"
}

# Create the Dash app
server = Flask(__name__)
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    external_scripts=[FONT_AWESOME],
    server=server,
)

app.title = "Bacterial Operon Finder for Functional Organization"


######################
# DOWNLOAD PDF ROUTE #
######################

# Folder used to host files
DOWNLOAD_DIR = os.path.join(os.getcwd(), "boffo_downloads")
if not os.path.exists(DOWNLOAD_DIR):
    os.mkdir(DOWNLOAD_DIR)

# Flask route
@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the download directory."""
    logger.info(path)
    return send_file(os.path.join(DOWNLOAD_DIR, path), as_attachment=True)

#####################
# SET UP APP LAYOUT #
#####################

def center_style(**kwargs):
    return {
        **kwargs,
        **dict(
            marginLeft="auto",
            marginRight="auto",
        )
    }


#######################
# PLOTTING PARAMETERS #
#######################

def parameter_form(
    name="Parameter Name",
    elem_id="element-id",
    options=[],
    value=True,
    inline=True,
    selector_type="radio",
    min_val=None,
    max_val=None,
    suffix="",
    step=1
):

    # Set the ID as a dict
    elem_dict = dict(elem_type="input", elem=elem_id)

    # If the element is the default, set up a radio selector
    if selector_type == "radio":

        selector = dbc.RadioItems(
            options=options,
            value=value,
            inline=inline,
            id=elem_dict
        )

    # If the element is a dropdown
    if selector_type == "dropdown":

        selector = dcc.Dropdown(
            options=options,
            value=value,
            id=elem_dict,
        )

    # If the element is a slider
    elif selector_type == "slider":

        # Mark the middle value
        mid_val = min_val + ((max_val - min_val) / 2.)

        selector = dcc.Slider(
            min=min_val,
            max=max_val,
            value=value,
            marks={x: f"{x}{suffix}" for x in [
                int(min_val) if min_val == int(min_val) else min_val, 
                int(mid_val) if mid_val == int(mid_val) else mid_val, 
                int(max_val) if max_val == int(max_val) else max_val
            ]},
            step=step,
            id=elem_dict,
        )

    # Format the entire card
    return dbc.FormGroup(
        [
            dbc.Col(
                [
                    dbc.Label(name),
                    selector,
                ],
                className="text-center"
            )
        ]
    )

param_list = [
    dict(
        name="Color Tree By Operon Presence",
        elem_id="color_by_operon",
        options=[
            dict(label="Yes", value=True),
            dict(label="No", value=False),
        ],
        value=True,
    ),
    dict(
        name="Show Genome Names",
        elem_id="genome_names",
        options=[
            dict(label="Yes", value=True),
            dict(label="No", value=False),
        ],
        value=True,
    ),
    dict(
        name="Show Genome IDs",
        elem_id="genome_ids",
        options=[
            dict(label="Yes", value=True),
            dict(label="No", value=False),
        ],
        value=False,
    ),
    dict(
        name="Show Gene Table",
        elem_id="show_genes",
        options=[
            dict(label="Yes", value=True),
            dict(label="No", value=False),
        ],
        value=True
    ),
    dict(
        name="Color Genes by Identity",
        elem_id="color_by_ident",
        options=[
            dict(label="Yes", value=True),
            dict(label="No", value=False),
        ],
        value=True
    ),
    dict(
        name="Minimum Alignment Identity",
        elem_id="min_iden",
        selector_type="slider",
        min_val=50,
        max_val=100,
        value=75,
        step=0.1,
        suffix="%"
    ),
    dict(
        name="Gene Table Width",
        elem_id="gene_table_width",
        selector_type="slider",
        min_val=0.01,
        max_val=1.0,
        value=0.5,
        step=0.01
    ),
    dict(
        name="Line Thickness",
        elem_id="line_width",
        selector_type="slider",
        min_val=0,
        max_val=5,
        value=1,
        step=0.1,
        suffix="px"
    ),
    dict(
        name="Figure Height",
        elem_id="figure_height",
        selector_type="slider",
        min_val=400,
        max_val=1200,
        value=600,
        suffix="px"
    ),
    dict(
        name="Figure Width",
        elem_id="figure_width",
        selector_type="slider",
        min_val=400,
        max_val=1200,
        value=800,
        suffix="px"
    ),
]

# Function to split up a list
def split_list(input_list, sublist_len):
    while len(input_list) > 0:
        yield input_list[:sublist_len]
        input_list = input_list[sublist_len:]


# Render a button using an icon
def button_icon(
    icon_text,
    icon_class_name,
    id_string,
    color="primary",
    outline=False,
    style=dict(
        marginLeft="10px"
    ),
    mouseover=None
):
    button_elem = dbc.Button(
        html.Span([
            icon_text,
            html.I(
                className=icon_class_name,
                style=dict(marginLeft="5px")
            )
        ]),
        id=id_string,
        color=color,
        outline=outline,
        style=style
    )

    if mouseover is None:
        return button_elem

    else:

        return html.Div([
            html.Div(
                button_elem,
                id=f"{id_string}-div"
            ),
            dbc.Tooltip(
                mouseover,
                target=f"{id_string}-div",
                placement="bottom"
            )
        ])


def custom_menu_collapse():
    # SHOW SETTINGS IN COLLAPSE
    return dbc.Collapse(
        dbc.Card(
            dbc.CardBody(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    parameter_form(
                                        **param
                                    )
                                    for param in param_sublist
                                ]
                            )
                            for param_sublist in split_list(param_list, 5)
                        ]
                    )
                ]
            )
        ),
        id='settings-collapse',
        is_open=False,
        style=dict(marginTop="10px")
    )


# Set up the layout of the app
app.layout = dbc.Container(
    [
        # Notification Toast
        dbc.Toast(
            "",
            id=toast_id,
            header="BOFFO",
            is_open=False,
            dismissable=True,
            icon="primary",
            style={
                "position": "fixed",
                "top": 92,
                "right": 10,
                "width": 350,
                "zIndex": 10,
            },
        )
        for toast_id in ["notification-toast", "download-notification-toast"]
    ] + [
        # THIS CARD IS THE HEADER ABOVE THE PLOT
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H4("Bacterial Operon Display")
                ),
                dbc.CardBody(
                    [
                        dbc.Label("Select Genes to Display"),
                        dcc.Dropdown(
                            options=[
                                dict(
                                    label="All Genes",
                                    value="ALL"
                                )
                            ] + [
                                    dict(label=v, value=v)
                                    for v in operon_list
                                ],
                            value="ALL",
                            id=dict(elem_type="input", elem="operon_context"),
                        ),
                        dbc.Form(
                            [
                                dbc.Spinner(
                                    dbc.FormGroup(
                                        dbc.Collapse(
                                            dbc.Button(
                                                button_text,
                                                id=button_id,
                                                color="primary",
                                                style=dict(
                                                    marginRight="10px",
                                                    marginTop="10px",
                                                    marginBottom="10px",
                                                ),
                                                href=""
                                            ),
                                            id=f"{button_id}-collapse",
                                            is_open=is_open,
                                        )
                                    )
                                )
                                for button_text, button_id, is_open in [
                                    ("Show Settings", "show-hide-settings", True),
                                    ("Make PDF", "make-pdf", True),
                                    ("Download PDF", "download-pdf", False),
                                    (html.Span(html.I(className="fas fa-redo")), "refresh-plot", True),
                                ]
                            ],
                            inline=True
                        )
                    ]
                )
            ],
            style=dict(marginTop="10px")
        ),
        # CARD WITH CUSTOMIZABLE PARAMETERS
        custom_menu_collapse(),
        # THE MAIN PLOT
        dbc.Row(
            html.Div(
                dcc.Graph(id="custom-plot"),
                style=center_style(marginTop="5px")
            )
        ),
        # THE LEGEND TEXT
        dbc.Row(
            html.Div(
                dcc.Markdown(id="plot-help-text"),
                style=center_style(
                    marginTop="5px",
                    marginLeft="20px",
                    marginRight="20px",
                    textAlign="justify",
                    textJustify="inter-word",
                )
            )
        )
    ]
)


####################
# SET UP CALLBACKS #
####################

# Decorate the callback functions with @app.callback as appropriate

# SHOW / HIDE SETTINGS
@app.callback(
    [
        Output("settings-collapse", "is_open"),
        Output("show-hide-settings", "children")
    ],
    [
        Input("show-hide-settings", "n_clicks")
    ],
    [
        State("settings-collapse", "is_open")
    ]
)
def show_hide_settings(n_clicks, is_open):
    if n_clicks is None:
        raise PreventUpdate
    elif is_open:
        return False, "Show Settings"
    else:
        return True, "Hide Settings"
    

# RENDER THE PLOT
@app.callback(
    [
        Output("custom-plot", "figure"),
        Output("plot-help-text", "children"),
        Output("notification-toast", "children"),
        Output("notification-toast", "is_open"),
    ],
    [
        Input({"elem_type": "input", "elem": ALL}, "value"),
        Input("refresh-plot", "n_clicks")
    ]
)
def render_plot_callback(
    input_list, _
):

    # Get the callback context in order to parse the inputs' IDs
    ctx = dash.callback_context

    # Parse the inputs to the callback
    input_values = {
        elem["id"]["elem"]: elem["value"]
        for elem in ctx.inputs_list[0]
    }

    # Try to render the plot, catching any errors
    try:
        fig = BOFFO_Plot(
            input_values
        ).fig

    # If there is an error
    except Exception as e:

        # Report the error to the user by opening the toast
        return report_error(str(e))

    # Try to render the help text, catching any errors
    try:
        text = render_text(
            input_values
        )

    # If there is an error
    except Exception as e:

        # Report the error to the user by opening the toast
        return report_error(str(e))

    # Otherwise, return the plot and don't open the toast
    return fig, text, "", False

# SAVE A PDF
@app.callback(
    [
        Output("make-pdf-collapse", "is_open"),
        Output("download-pdf-collapse", "is_open"),
        Output("download-pdf", "href"),
        Output("download-pdf", "external_link"),
        Output("download-notification-toast", "children"),
        Output("download-notification-toast", "is_open"),
    ],
    [
        Input("download-pdf", "n_clicks"),
        Input("make-pdf", "n_clicks"),
    ],
    [
        State({"elem_type": "input", "elem": ALL}, "value"),
    ]
)
def download_plot_callback(
    download_clicks, make_pdf_clicks, input_list
):

    logger.info("download_plot_callback")

    # Get the callback context in order to parse the inputs' IDs
    ctx = dash.callback_context

    # Figure out which button triggered the callback
    trigger_elem = ctx.triggered[0]['prop_id']
    
    # If there was no trigger
    if trigger_elem == ".":

        # Take no action
        raise PreventUpdate

    # If the Download PDF button was pressed
    elif trigger_elem == "download-pdf.n_clicks":

        # Hide the download button and reveal the Make PDF button
        return True, False, "", False, "", False

    # If the Make PDF button was pressed
    elif trigger_elem == "make-pdf.n_clicks":

        # Parse the inputs to the callback
        input_values = {
            elem["id"]["elem"]: elem["value"]
            for elem in ctx.states_list[0]
        }

        # Try to render the plot, catching any errors
        try:
            fig = BOFFO_Plot(
                input_values
            ).fig

        # If there is an error
        except Exception as e:

            # Report the error to the user by opening the toast
            return True, False, "", False, str(e), True

        # Otherwise, the figure was made

        # Make a random string for the download
        pdf_fn = f"{str(uuid4())[:8]}.pdf"

        logger.info(pdf_fn)

        # Save it to a PDF
        try:

            # Set up the path to the image
            img_path = os.path.join(DOWNLOAD_DIR, pdf_fn)

            # Write the image
            fig.write_image(img_path)

            # Make sure that it exists
            assert os.path.exists(img_path)

        # If there is an error
        except Exception as e:

            # Report the error to the user by opening the toast
            return True, False, "", False, str(e), True

        # Otherwise, populate the link to the file via the flask download route
        return False, True, f"/download/{pdf_fn}", True, "", False


def line_to_parent(node, line_width):
    """Draw a line from a node to its parent."""

    # Get the X-Y of the parent
    parent_x = tree_layout.nodes[node["parent"]]["x"]
    parent_y = tree_layout.nodes[node["parent"]]["y"]

    # Format the X and Y coordinates
    x_pos = [
        node["x"],
        parent_x,
        parent_x
    ]
    y_pos = [
        node["y"],
        node["y"],
        parent_y
    ]

    # Return a line
    return go.Scatter(
        x=x_pos,
        y=y_pos,
        mode="lines",
        showlegend=False,
        marker_color="black",
        line_width=float(line_width)
    )


@lru_cache(maxsize=128)
def get_genome_operon_labels(operon_context):
    """For each genome, label it based on whether the operon is found."""

    # Make an empty dict
    output = dict()

    # If "Show All Genes" was selected
    if operon_context == "ALL":

        # Then fill in a "None" value for each genome
        return output

    # Iterate over each genome which has the operon
    for genome_id in operons.query(
        f"operon_context == '{operon_context}'"
    )[
        "genome_id"
    ].values:

        # Mark it as being present
        output[genome_id] = "Complete Operon"

    # Get the unique list of genes present in this operon
    genes_in_operon = operons.query(
        f"operon_context == '{operon_context}'"
    )[
        "gene_name"
    ].unique()

    # Iterate over every genome which contains any of those genes
    for genome_id in operons.loc[
        operons["gene_name"].isin(genes_in_operon)
    ][
        "genome_id"
    ].unique():

        # If the genome is not already labelled
        if output.get(genome_id) is None:

            # Then label the genome "Partial Operon"
            output[genome_id] = "Partial Operon"

    # NOTE: Genomes with no genes detected have empty keys

    # Return the dict
    return output

class BOFFO_Plot:
    """Set up the main figure."""
    
    def __init__(
        self,
        params
    ):
        """Render the main figure."""

        # Attach the plot parameters to the object
        self.params = params

        # Initialize the figure
        self.fig = self.setup_figure()

        # Set up objects to hold the coordinates
        self.x_pos = []      # Store the x coordinates for each index position
        self.y_pos = []      # Store the y coordinates for each index position
        self.names = []      # Store the name to display at each index position
        self.ids = []        # Store the unique genome ID for each index position
        self.name_dict = {}  # Map the unique genome ID to the selected genome name

        # Set up the position of all genomes
        self.plot_genome_tree()

        # If the user elected to color the genomes by operon presence
        if self.params["color_by_operon"]:

            self.color_points_by_operon()

        # If the user elected to show the gene presence-absence table
        if self.params["show_genes"]:

            self.add_gene_heatmap(
                color_by_ident=self.params["color_by_ident"],
                min_iden=self.params["min_iden"]
            )

        # Finally, add the genome labels
        self.add_genome_labels()


    def setup_figure(self):
        
        # Set up a plotly figure
        fig = go.Figure()

        # Update the layout of the figure
        fig.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            yaxis_showticklabels=True,
            xaxis_title="Branch Length",
            height=int(self.params["figure_height"]),
            width=int(self.params["figure_width"]),
            legend=dict(
                yanchor="bottom",
                y=1,
                xanchor="left",
                x=0,
                orientation="v",
            ),
            showlegend=True,
            coloraxis_showscale=False,
        )

        return fig

    def plot_genome_tree(self):

        # Get the description of each genome
        # as either "Complete Operon", "Partial Operon", or None
        genome_operon_labels = get_genome_operon_labels(
            self.params["operon_context"]
        )

        # Make a dict with the positions of each group of genomes
        self.operon_xy_dict = {
            v: dict(
                x=[],
                y=[]
            )
            for v in list(set(list(genome_operon_labels.values())))
        }

        # Iterate over every node
        for node_name, node in tree_layout.nodes.items():

            # Record the X and Y position for the node (by index position)
            self.x_pos.append(node["x"])
            self.y_pos.append(node["y"])

            # Record the unique genome ID for this node (by index position)
            self.ids.append(node_name)

            # If there is a parent node
            if node["parent"] is not None:

                # Draw a line to the parent
                self.fig.add_trace(
                    line_to_parent(
                        node,
                        self.params["line_width"]
                    )
                )

            # If this is a terminal node
            if node["terminal"]:

                # If the genome ID -> name mapping failed
                if genome_name_dict.get(node_name) is None:

                    # Then just add the ID
                    self.names.append(node_name)

                # Otherwise, if the user elected to show both the name and ID
                elif self.params["genome_names"] and self.params["genome_ids"]:

                    self.names.append(
                        f"{genome_name_dict[node_name]} ({node_name})"
                    )

                # If the user just elected to show the name
                elif self.params["genome_names"]:
                    self.names.append(genome_name_dict[node_name])

                # If the user just elected to show the name
                elif self.params["genome_ids"]:

                    # Just add the ID
                    self.names.append(node_name)

                # Otherwise
                else:

                    # Add an empty string
                    self.names.append("")

                # Add the formatted name to the dict keyed on unique ID
                self.name_dict[node_name] = self.names[-1]
                
                # If this genome contains the operon
                if genome_operon_labels.get(node_name) is not None:

                    # Add the point to the operon_xy_dict dict
                    for v in ["x", "y"]:
                        self.operon_xy_dict[
                            genome_operon_labels.get(node_name)
                        ][v].append(node[v])

            # Otherwise, if it is internal
            else:

                # Add a blank string
                self.names.append("")

        # Add a point for each node in the tree
        self.fig.add_trace(
            go.Scatter(
                x=self.x_pos,
                y=self.y_pos,
                mode="markers",
                showlegend=False,
                marker_color="black",
                cliponaxis=False,
                hoverinfo="text",
                hovertext=self.names
            )
        )

    def color_points_by_operon(self):

        # For each type of operon presence
        for label, vals in self.operon_xy_dict.items():

            # Add the points to the plot
            self.fig.add_trace(
                go.Scatter(
                    x=vals["x"],
                    y=vals["y"],
                    mode="markers",
                    showlegend=True,
                    hoverinfo="skip",
                    name=label
                )
            )

    def add_genome_labels(self):

        # Add genome labels as yticklabels
        
        # Make a list of the y-positions and genome names
        ytick_pos = []
        ytick_label = []

        # Iterate over each of the names and y-coordinates in the tree
        for name, y in zip(self.names, self.y_pos):

            # If the node is not internal
            if len(name) > 0:

                # Record this particular position and name
                ytick_pos.append(y)
                ytick_label.append(name)

        # Update the figure
        self.fig.update_layout(
            yaxis = dict(
                side="right",
                tickmode = 'array',
                tickvals = ytick_pos,
                ticktext = ytick_label
            )
        )

    def cluster_genes_by_presence(self, metric="braycurtis", method="average"):
        """Perform linkage clusering on the genes detected."""

        wide_df = operons.pivot_table(
            index="gene_name",
            columns="genome_id",
            values="pct_iden",
            aggfunc=max
        ).fillna(
            0
        )

        return [
            wide_df.index.values[i]
            for i in leaves_list(
                linkage(
                    wide_df,
                    metric=metric,
                    method=method
                )
            )
        ]

    def add_gene_heatmap(self, color_by_ident=True, min_iden=0):

        # Get the operon context to display
        operon_context = self.params["operon_context"]

        # If "All Genes" was selected
        if operon_context == "ALL":

            # Show all genes
            operon_df = operons

            # Order the genes in the operon by hierarchical clustering
            genes_in_operon = self.cluster_genes_by_presence()
        
        # Otherwise, if a specific operon was selected
        else:
            
            # Get the list of genes which are in this operon
            genes_in_operon = [
                g[:-4]  # Remove the (+) or (-) from the end of the gene name
                for g in operon_context.split(" :: ")
            ]

            # Get the table for alignments to these genes
            operon_df = operons.loc[
                operons["gene_name"].isin(genes_in_operon)
            ]

        # Filter down the data to display
        operon_df = operon_df.drop(
            columns=[
                "aligned_sequence",
                "translated_sequence",
                "genome_context",
                "operon_size",
                "alignment_length",
                "gene_start",
                "gene_end",
                "gene_len",
            ]
        ).loc[
            operon_df["pct_iden"].apply(float) >= min_iden
        ]

        assert operon_df.shape[0] > 0, "No alignments pass the minimum identity filter"

        # Get the y position of every genome by id
        genome_y = dict(zip(self.ids, self.y_pos))

        # Set up the lists of points which will be plotted
        genes_x = []
        genes_y = []
        genes_hovertext = []
        pct_iden = []

        # Iterate over every genome in the table
        for genome_id, genome_df in operon_df.groupby("genome_id"):

            # Iterate over every gene which was aligned
            for gene_id, df in genome_df.groupby("gene_name"):

                # Format the text to plot
                genes_hovertext.append(self.format_hovertext(df))

                # Add the X coordinate for the gene
                genes_x.append(genes_in_operon.index(gene_id))

                # Add the Y coordinate for the genome
                genes_y.append(genome_y[genome_id])

                # Add the percent identity for the alignment
                pct_iden.append(df["pct_iden"].max())

        # Add the points to the plot
        self.fig.add_trace(
            go.Scatter(
                x=genes_x,
                y=genes_y,
                mode="markers",
                showlegend=False,
                marker=dict(
                    color=pct_iden if color_by_ident else "black",
                    colorscale="Bluered",
                ),
                hoverinfo="text",
                hovertext=genes_hovertext,
                xaxis="x2",
            )
        )

        # Left-align the hovertext
        # Add labels for the genes above the heatmap
        self.fig.update_layout(
            hoverlabel_align="left",
            xaxis2=dict(
                side="top",
                tickmode = 'array',
                tickvals = list(range(len(genes_in_operon))),
                ticktext = genes_in_operon,
                domain = [1 - self.params["gene_table_width"], 1],
            ),
            xaxis=dict(
                domain = [0, 1 - self.params["gene_table_width"]],
            ),
            yaxis=dict(
                anchor="x2"
            )
        )

    def format_hovertext(self, df):
        """Given a list of BOFFO alignments, format the aggregate hovertext."""

        # Get the genome name
        genome_name = self.name_dict[df["genome_id"].values[0]]

        # Get the gene name
        gene_name = df["gene_name"].values[0]
        
        # Format the alignment information for each alignment
        alignment_text = "<br>".join(
            df.apply(
                self.format_hovertext_row,
                axis=1
            )
        )

        # Format the entire string
        return f"{genome_name} ::: {gene_name}<br><br>{alignment_text}"

    def format_hovertext_row(self, r):
        """Format the hovertext for a single alignment."""
        cstart = r['contig_start']
        cend = r['contig_end']
        cname  = r["contig_name"]
        cov = round(r["gene_cov"], 1)
        pcti = round(r["pct_iden"], 1)
        return f"{cname} ({cstart:,} - {cend:,}): {cov}% gene length, {pcti}% amino acid identity"



def render_text(input_values):
    """Render the help text."""

    text = """
**Comparison of global genomic similarity with the presence of specific genes.**

The dendrogram on the left side of the plot represents a neighbor-joining tree 
constructed by [mashtree](https://github.com/lskatz/mashtree), which compares the
total nucleotide content across a collection of genomes. This analysis should not
necessarily be interpreted as representing evolutionary history because it does not
account for any confounding influence of inaccurate genome assemblies or horizontal
transfer of genetic material. The aggregate horizontal distance of each leaf
represents the proportion of genetic content which differs between any pair of genomes.
"""

    if input_values["color_by_operon"] and input_values["operon_context"] != "ALL":
        text = f"""{text}
Each genome has been annotated by color to indicate whether it contains
the complete operon (including all genes in the indicated relative orientation)
or a partial set of genes in the operon.
"""

    if input_values["show_genes"]:

        text = f"""
{text}
The points on the right-hand side of the plot indicate whether each genome
contains any of the selected genes, with a marker indicating presence.
Moving your cursor over each marker will display the coordinates of
each alignment for the indicated gene along with the coverage and percent identity
of the alignment.
        """

    if input_values["color_by_ident"]:
        text = f"""
{text}
The color of each marker reflects the percent identity of the alignment, with
lower values of amino acid similarity shown in red and higher values in blue.
"""

    text = f"""
{text}

Use the drop-down menu at the top of the page to select the operon to display.
Use the display settings menu below to customize the plot further. The final
image may be downloaded by clicking the "Make PDF" button, followed by "Download PDF".

The user may expand and zoom in on different parts of the dendrogram or gene
presence matrix by clicking and dragging on the plot directly. The plot position
can be reset by double-clicking on the plot. Changes to the plot made by zooming
in this way will not be reflected with the "Make PDF" button, but they can be
saved in SVG format by clicking on the camera icon within the plot area.
"""
    return text


def report_error(error_text):
    return go.Figure(), "", error_text, True


##############
# SET UP APP #
##############

# Run the app
if __name__ == '__main__':

    app.run_server(
        host='0.0.0.0',
        port=8080,
        debug=False
    )
