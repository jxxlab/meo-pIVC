# ============================================================
# Figure 3d/e/f RNA velocity plotting — condensed reusable code
# Summarized from:
# demo-scVelomoe-moe-sc-E9-D6-RD4-C1-C12-185873.html
#
# Logic:
#   1. Read one lineage/subset h5ad file containing spliced/unspliced layers.
#   2. Run scVelo preprocessing and velocity analysis.
#   3. Keep the same CX1 color mapping as the original Figure 3 palette.
#   4. Split the same dataset into:
#        Vivo   = in vivo
#        Vitro1 = meo-pIVC1
#        Vitro2 = meo-pIVC2
#   5. For EACH group, recalculate neighbors, moments, velocity and velocity graph.
#
# Figure 3d/e/f use the same workflow.
# For the other panels, only replace input_h5ad and output_prefix.
# ============================================================


# ============================================================
# 0. Packages and global settings
# ============================================================

import os
import numpy as np
import pandas as pd
import scvelo as scv
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PathCollection

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")

scv.set_figure_params(
    dpi=100,
    fontsize=16,
    dpi_save=100,
    figsize=(3, 3),
    frameon=False
)

# Editable text in PDF/SVG
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


# ============================================================
# 1. Original CX1 colors
# ============================================================

x_colors = [
    "#FCFF00", "#BD84B0", "#456722", "#FFAAFF", "#BF9D88", "#4A6798",
    "#FACB12", "#0F4A9C", "#FFD731", "#685ED8", "#E3CB3A", "#683ED8",
    "#2FDA00", "#A64D7E", "#AAFFAA", "#9D0049", "#5ADBE4", "#772600",
    "#89C1F5", "#BD3400", "#3F84AA", "#DDAA22", "#5581CA", "#DEA93C",
    "#2F4A60", "#F79083", "#127D4C", "#FF7F9C", "#2F9A00", "#B51D8D",
    "#549E79", "#E85639", "#00BFC4", "#FF5C00", "#005579", "#FBBE92",
    "#532C5A", "#E3CB9A", "#532C8A", "#DABE99", "#7C2A47", "#F9DFE6",
    "#635547", "#C594BF", "#D5E839", "#683EA8", "#C19F70"
]


# ============================================================
# 2. Build the CX1 palette
# ============================================================

def prepare_cx1_palette(adata):

    # CX1 -> string
    adata.obs["CX1"] = adata.obs["CX1"].astype(str)

    # Sort only clusters that really exist in this object
    cx1_categories = sorted(
        adata.obs["CX1"].unique(),
        key=lambda z: int(z)
    )

    adata.obs["CX1"] = pd.Categorical(
        adata.obs["CX1"],
        categories=cx1_categories,
        ordered=True
    )

    # Cluster number corresponds to the same position in x_colors
    full_palette = {
        str(i + 1): color.upper()
        for i, color in enumerate(x_colors)
    }

    x_palette = {
        cluster: full_palette[cluster]
        for cluster in cx1_categories
    }

    return x_palette


# ============================================================
# 3. scVelo analysis for ONE culture group
# ============================================================
# Important:
# The original notebook recalculated the neighborhood graph,
# moments, velocity and velocity graph separately for Vivo,
# Vitro1 and Vitro2.

def run_velocity_for_group(adata_all, group_name):

    adata = adata_all[
        adata_all.obs["sample2"].isin([group_name])
    ].copy()

    print(group_name, adata.shape)

    # RNA velocity preprocessing
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=20,
        n_top_genes=2000
    )

    # Recalculate the group's own neighborhood graph
    scv.pp.neighbors(
        adata,
        n_pcs=30,
        n_neighbors=30
    )

    # Recalculate moments
    scv.pp.moments(
        adata,
        n_pcs=30,
        n_neighbors=30
    )

    # Steady-state velocity model used in the notebook
    scv.tl.velocity(adata)

    # Velocity graph
    scv.tl.velocity_graph(
        adata,
        n_jobs=10
    )

    return adata


# ============================================================
# 4. Common plotting function
# ============================================================

def plot_velocity_group(
    adata,
    x_palette,
    group_name,
    output_prefix
):

    fig, ax = plt.subplots(figsize=(4, 4))

    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="CX1",
        palette=x_palette,

        size=10,
        alpha=0.5,

        density=0.6,
        arrow_size=1,
        linewidth=1.0,
        smooth=0.8,
        max_length=6,

        legend_loc="on data",
        legend_fontsize=6,
        legend_fontweight="bold",

        title="",
        ax=ax,
        show=False
    )

    # --------------------------------------------------------
    # Keep cluster labels editable
    # --------------------------------------------------------

    for txt in ax.texts:
        txt.set_path_effects([])
        txt.set_rasterized(False)
        txt.set_fontsize(6)
        txt.set_fontweight("bold")
        txt.set_color("black")

    # --------------------------------------------------------
    # Fix possible NaN / Inf linewidth values
    # --------------------------------------------------------

    for i, artist in enumerate(ax.collections):

        if hasattr(artist, "get_linewidths"):

            lw = np.asarray(
                artist.get_linewidths(),
                dtype=float
            )

            if lw.size > 0:

                n_bad = np.sum(~np.isfinite(lw))

                if n_bad > 0:

                    print(
                        f"Fix collection {i}: "
                        f"{type(artist).__name__}, "
                        f"{n_bad} non-finite linewidths"
                    )

                    lw = np.nan_to_num(
                        lw,
                        nan=0.0,
                        posinf=0.0,
                        neginf=0.0
                    )

                    artist.set_linewidths(lw)

    # --------------------------------------------------------
    # Rasterize cell points only.
    # Streamlines and text remain vector/editable.
    # --------------------------------------------------------

    for artist in ax.collections:
        if isinstance(artist, PathCollection):
            artist.set_rasterized(True)

    # --------------------------------------------------------
    # Save editable PDF
    # --------------------------------------------------------

    output_file = (
        f"{output_prefix}_{group_name}_editable.pdf"
    )

    fig.savefig(
        output_file,
        format="pdf",
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )

    plt.show()
    plt.close(fig)

    print("Saved:", output_file)


# ============================================================
# 5. One complete Figure 3 panel
# ============================================================
# One input h5ad produces the three plots:
#
#   Vivo   -> in vivo
#   Vitro1 -> meo-pIVC1
#   Vitro2 -> meo-pIVC2
#
# d/e/f can all call this same function.

def make_velocity_panel(
    input_h5ad,
    output_prefix
):

    # --------------------------------------------------------
    # Read lineage/subset h5ad
    # --------------------------------------------------------

    adata = scv.read(input_h5ad)

    print(adata)

    # --------------------------------------------------------
    # The notebook first performed preprocessing on the
    # complete lineage/subset dataset.
    # --------------------------------------------------------

    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=20,
        n_top_genes=2000
    )

    scv.pp.moments(
        adata,
        n_pcs=30,
        n_neighbors=30
    )

    scv.tl.velocity(adata)

    scv.tl.velocity_graph(
        adata,
        n_jobs=10
    )

    # --------------------------------------------------------
    # Prepare colors once from the whole lineage object
    # --------------------------------------------------------

    x_palette = prepare_cx1_palette(adata)

    print(x_palette)

    # Preserve the complete object before splitting
    adata_all = adata.copy()

    # --------------------------------------------------------
    # Generate the three matched plots
    # --------------------------------------------------------

    for group_name in [
        "Vivo",
        "Vitro1",
        "Vitro2"
    ]:

        adata_group = run_velocity_for_group(
            adata_all,
            group_name
        )

        plot_velocity_group(
            adata_group,
            x_palette,
            group_name,
            output_prefix
        )


#
# /sdc/jxx2/moe/moe-sc-E9-D6-RD4-C1-C12-185873.h5ad
#
# This object contained the C1-C12 subset.

make_velocity_panel(
    input_h5ad=(
        "/sdc/jxx2/moe/"
        "moe-sc-E9-D6-RD4-C1-C12-185873.h5ad"
    ),
    output_prefix=(
        "velocity_moe-sc-E9-D6-RD4-C1-C12-185873"
    )
)

