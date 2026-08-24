#!/usr/bin/env python3
"""
============================================================================
 MD / Desmond-SID publication figure generator
============================================================================
Generates all six standard Simulation-Interaction-Diagram panels from the
raw .dat exports, in a consistent publication style:

    PL_RMSD.dat          -> PL_RMSD.png/.pdf     (protein-ligand RMSD vs time)
    P_RMSF.dat           -> P_RMSF.png/.pdf       (protein RMSF vs residue)
    L_RMSF.dat           -> L_RMSF.png/.pdf       (ligand per-atom RMSF)
    L-Properties.dat     -> L_Properties.png/.pdf (6 ligand descriptors vs time)
    L_Torsions.dat       -> L_Torsions.png/.pdf   (rotatable-bond dial plots)
    PL-Contacts_*.dat    -> PL_Contacts.png/.pdf  (interaction-fraction bars)

----------------------------------------------------------------------------
HOW TO REUSE ON A NEW RESULT SET
----------------------------------------------------------------------------
Put a result set's .dat files (same file names) in one folder, then:

    python3 plot_md_sid.py --indir  path/to/run2 \
                           --outdir path/to/run2_figures \
                           --time   200

  --indir   folder containing the .dat files      (default: current folder)
  --outdir  folder to write the figures into       (default: same as --indir)
  --time    total simulation length in ns          (default: 200)

The number of frames is detected automatically from PL_RMSD.dat /
L-Properties.dat, so runs of different length work without edits.
Only --time must be set correctly for each run (it is not stored in the
.dat files). Missing or empty inputs are skipped with a message.
============================================================================
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from collections import Counter

# ----------------------------------------------------------------------
# STYLE (shared by every figure)
# ----------------------------------------------------------------------
UNIT = "\u00c5"                       # Angstrom symbol
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,

    "axes.labelsize":16,
    "axes.titlesize":18,
    "axes.titleweight":"bold",

    "xtick.labelsize":13,
    "ytick.labelsize":13,

    "legend.fontsize":12,

    "axes.linewidth":1.4,

    "xtick.direction":"in",
    "ytick.direction":"in",

    "xtick.major.size":6,
    "ytick.major.size":6,

    "xtick.minor.size":3,
    "ytick.minor.size":3,

    "savefig.dpi":600,
    "figure.dpi":150
})
C_PROT  = "#1f4e79"   # blue   - protein
C_LIG   = "#c0392b"   # red    - ligand
C_LINE  = "#5b7488"   # slate  - single traces
C_HELIX = "#e8b4b4"   # soft red  - alpha-helix band
C_STRAND= "#b4cfe0"   # soft blue - beta-strand band
C_CONT  = "#1a7a34"   # green  - ligand-contact bars

# ----------------------------------------------------------------------
# OPTIONAL: protein secondary-structure bands for the P_RMSF panel.
# Not present in P_RMSF.dat (comes from the Desmond SSE export or DSSP).
# Fill as (start_index, end_index) pairs in RESIDUE-INDEX units, or leave
# empty to omit. These are per-run, so edit here when you have them.
# ----------------------------------------------------------------------
HELIX_RANGES  = []    # e.g. [(12, 24), (58, 71)]
STRAND_RANGES = []    # e.g. [(3, 9), (40, 46)]


# ======================================================================
# helpers
# ======================================================================
def load_numeric(path):
    """Load whitespace-separated numeric rows, skipping any header line
    (handles '#' at line start OR mid-line such as 'Frame #')."""
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return np.array(rows)


def detect_nframes(indir):
    """Number of trajectory frames, from PL_RMSD.dat or L-Properties.dat."""
    for fn in ("PL_RMSD.dat", "L-Properties.dat"):
        p = os.path.join(indir, fn)
        if os.path.exists(p):
            d = load_numeric(p)
            if d.size:
                return d.shape[0]
    return None


def _save(fig, outdir, name):
    png = os.path.join(outdir, name + ".png")
    pdf = os.path.join(outdir, name + ".pdf")
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.png / .pdf")


# ======================================================================
# 1. PROTEIN-LIGAND RMSD  (dual axis vs time)
# ======================================================================
def plot_rmsd(indir, outdir, total_time,
              protein_col="Prot_CA", ligand_col="Lig_wrt_Protein"):
    path = os.path.join(indir, "PL_RMSD.dat")
    if not os.path.exists(path):
        print("  [skip] PL_RMSD.dat not found"); return
    names = ["frame", "Prot_CA", "Prot_Backbone", "Prot_Sidechain",
             "Prot_All_Heavy", "Lig_wrt_Protein", "Lig_wrt_Ligand"]
    col = {n: i for i, n in enumerate(names)}
    d = load_numeric(path)
    t = d[:, 0] / d[-1, 0] * total_time
    prot, lig = d[:, col[protein_col]], d[:, col[ligand_col]]

    fig, ax1 = plt.subplots(figsize=(8,5))
    ax1.plot(t, prot, color=C_PROT, lw=2.0, alpha=0.9)
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel(f"Protein RMSD ({UNIT})", color=C_PROT)
    ax1.tick_params(axis="y", colors=C_PROT); ax1.spines["left"].set_color(C_PROT)

    ax2 = ax1.twinx()
    ax2.plot(t, lig, color=C_LIG, lw=2.0, alpha=0.9)
    ax2.set_ylabel(f"Ligand RMSD ({UNIT})", color=C_LIG)
    ax2.tick_params(axis="y", colors=C_LIG)
    ax2.spines["right"].set_color(C_LIG); ax2.spines["left"].set_visible(False)

    ax1.set_xlim(0, total_time)
    ax1.xaxis.set_major_locator(MultipleLocator(max(1, round(total_time/4))))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ymax = max(prot.max(), lig.max()) * 1.15
    ax1.set_ylim(0, ymax); ax2.set_ylim(0, ymax)
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax1.legend([Line2D([0], [0], color=C_PROT, lw=1.6),
                Line2D([0], [0], color=C_LIG, lw=1.6)],
               [f"Protein ({protein_col.replace('Prot_','')})",
                f"Ligand ({ligand_col.replace('Lig_wrt_','wrt ')})"],
               loc="upper left", frameon=False, fontsize=13)
    eq = t >= 0.25 * total_time
    ax1.text(0.985, 0.04,
             f"Protein: {prot[eq].mean():.2f} \u00b1 {prot[eq].std():.2f} {UNIT}\n"
             f"Ligand:  {lig[eq].mean():.2f} \u00b1 {lig[eq].std():.2f} {UNIT}",
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=10,
             color="#333", bbox=dict(boxstyle="round,pad=0.4", fc="white",
                                     ec="#ccc", alpha=0.9))
    _save(fig, outdir, "PL_RMSD")


# ======================================================================
# 2. PROTEIN RMSF  (vs residue index, SID style)
# ======================================================================
def plot_protein_rmsf(indir, outdir, rmsf_col="CA"):
    path = os.path.join(indir, "P_RMSF.dat")
    if not os.path.exists(path):
        print("  [skip] P_RMSF.dat not found"); return
    colidx = {"CA": 4, "Backbone": 5, "Sidechain": 6, "All_Heavy": 7, "B-factor": 8}
    res_index, contact, vals = [], [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) < 5:
                continue
            res_index.append(int(p[0]))
            contact.append(p[3] == "Yes")
            vals.append(float(p[colidx[rmsf_col]]))
    res_index = np.array(res_index); contact = np.array(contact); vals = np.array(vals)

    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    ymax = vals.max() * 1.10
    for a, b in HELIX_RANGES:
        ax.axvspan(a - .5, b + .5, color=C_HELIX, alpha=.55, lw=0, zorder=0)
    for a, b in STRAND_RANGES:
        ax.axvspan(a - .5, b + .5, color=C_STRAND, alpha=.55, lw=0, zorder=0)
    bar_top = 0.16 * ymax
    for r in res_index[contact]:
        ax.plot([r, r], [0, bar_top], color=C_CONT, lw=2,
                solid_capstyle="butt", zorder=1)
    ax.plot(res_index, vals, color=C_LINE, lw=2, zorder=3)

    ax.set_title("Protein RMSF", fontsize=14, fontweight="bold", pad=18)
    ax.set_xlabel("Residue Index"); ax.set_ylabel(f"RMSF ({UNIT})")
    ax.set_xlim(res_index.min(), res_index.max()); ax.set_ylim(0, ymax)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(0.6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    handles = [
        Line2D(
            [0], [0],
            color=C_LINE,
            lw=2,
            label="Cα" if rmsf_col == "CA" else rmsf_col
        )
    ]

    if HELIX_RANGES:
        handles.append(
            Patch(facecolor=C_HELIX,
                  alpha=0.55,
                  label="α-helix")
        )

    if STRAND_RANGES:
        handles.append(
            Patch(facecolor=C_STRAND,
                  alpha=0.55,
                  label="β-strand")
        )

    if contact.any():
        handles.append(
            Line2D(
                [0], [0],
                color=C_CONT,
                lw=2,
                label="Ligand Contact"
            )
        )

    ax.legend(
        handles=handles,
        loc="upper right",
        frameon=False,
        fontsize=12
    )

    fig.tight_layout()

    _save(fig, outdir, "P_RMSF")
# ======================================================================
# 3. LIGAND RMSF (per atom)
# ======================================================================

def plot_ligand_rmsf(indir, outdir):

    path = os.path.join(indir, "L_RMSF.dat")

    if not os.path.exists(path):
        print("  [skip] L_RMSF.dat not found")
        return

    d = load_numeric(path)

    atom = d[:, 0].astype(int)
    fit_protein = d[:, -2]
    fit_ligand = d[:, -1]

    fig, ax = plt.subplots(figsize=(8, 5))

    # ---------------------------
    # Plot curves
    # ---------------------------
    ax.plot(
        atom,
        fit_protein,
        color=C_PROT,
        linewidth=2.2,
        marker="o",
        markersize=5,
        label="Fit on Protein"
    )

    ax.plot(
        atom,
        fit_ligand,
        color=C_LIG,
        linewidth=2.2,
        marker="s",
        markersize=5,
        label="Fit on Ligand"
    )

    # ---------------------------
    # Titles and labels
    # ---------------------------
    ax.set_title(
        "Ligand RMSF",
        fontsize=18,
        fontweight="bold",
        pad=16
    )

    ax.set_xlabel(
        "Atom Index",
        fontsize=16
    )

    ax.set_ylabel(
        f"RMSF ({UNIT})",
        fontsize=16
    )

    # ---------------------------
    # Axis limits
    # ---------------------------
    ax.set_xlim(atom.min() - 0.5, atom.max() + 0.5)
    ax.set_ylim(0, max(fit_protein.max(), fit_ligand.max()) * 1.10)

    # ---------------------------
    # Tick formatting
    # ---------------------------
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=13,
        length=6,
        width=1.2
    )

    ax.tick_params(
        axis="both",
        which="minor",
        length=3,
        width=1.0
    )

    # ---------------------------
    # Legend
    # ---------------------------
    ax.legend(
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=12,
        handlelength=2.5,
        columnspacing=1.5
    )

    # ---------------------------
    # Final layout
    # ---------------------------
    fig.tight_layout()

    _save(fig, outdir, "L_RMSF")

# ======================================================================
# 4. LIGAND PROPERTIES (6 time series)
# ======================================================================

def plot_ligand_properties(indir, outdir, total_time):

    path = os.path.join(indir, "L-Properties.dat")

    if not os.path.exists(path):
        print("  [skip] L-Properties.dat not found")
        return

    d = load_numeric(path)

    t = d[:, 0] / d[-1, 0] * total_time

    properties = [
        ("RMSD",    d[:,1], f"RMSD ({UNIT})"),
        ("Radius of Gyration",    d[:,2], f"Radius of Gyration ({UNIT})"),
        ("intraHB", d[:,3], "Intramolecular H-bonds"),
        ("MolSA",   d[:,4], f"MolSA ({UNIT}²)"),
        ("SASA",    d[:,5], f"SASA ({UNIT}²)"),
        ("PSA",     d[:,6], f"PSA ({UNIT}²)")
    ]

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(10,8),
        sharex=True
    )

    fig.suptitle(
        "Ligand Properties",
        fontsize=20,
        fontweight="bold",
        y=0.98
    )

    for ax, (title, values, ylabel) in zip(axes.ravel(), properties):

        ax.plot(
            t,
            values,
            color=C_LINE,
            linewidth=1.8
        )

        ax.set_ylabel(
            ylabel,
            fontsize=13
        )

        ax.set_xlim(0, total_time)

        ax.set_title(
            title,
            fontsize=14,
            fontweight="bold",
            pad=8
        )

        ax.xaxis.set_major_locator(
            MultipleLocator(max(1, round(total_time/4)))
        )

        ax.xaxis.set_minor_locator(
            AutoMinorLocator(5)
        )

        ax.yaxis.set_minor_locator(
            AutoMinorLocator(2)
        )

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=11,
            length=5,
            width=1.2
        )

        ax.tick_params(
            axis="both",
            which="minor",
            length=3,
            width=1.0
        )

        ax.margins(y=0.10)

    # Bottom x-axis labels only
    for ax in axes[-1, :]:
        ax.set_xlabel(
            "Time (ns)",
            fontsize=13
        )

    fig.tight_layout(rect=[0,0,1,0.96])

    _save(fig, outdir, "L_Properties")
# ======================================================================
# 5. LIGAND TORSION PROFILES (Polar Histograms)
# ======================================================================

def plot_ligand_torsions(indir, outdir, nbins=36):

    path = os.path.join(indir, "L_Torsions.dat")

    if not os.path.exists(path):
        print("  [skip] L_Torsions.dat not found")
        return

    d = load_numeric(path)

    angles = d[:,1:]

    n = angles.shape[1]

    ncol = 3
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(10, 3.5*nrow),
        subplot_kw={"projection":"polar"}
    )

    fig.suptitle(
        "Ligand Torsion Profiles",
        fontsize=20,
        fontweight="bold",
        y=0.98
    )

    axes = np.atleast_1d(axes).ravel()

    bins = np.linspace(-np.pi, np.pi, nbins + 1)

    width = np.diff(bins)

    centers = bins[:-1] + width/2

    for i in range(n):

        ax = axes[i]

        counts, _ = np.histogram(
            np.deg2rad(angles[:,i]),
            bins=bins
        )

        ax.bar(
            centers,
            counts,
            width=width,
            color=C_LINE,
            alpha=0.90,
            edgecolor="white",
            linewidth=0.5
        )

        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)

        ax.set_thetagrids(
            [0,90,180,270],
            labels=["0","90","180","−90"],
            fontsize=10
        )

        ax.set_yticklabels([])

        ax.set_title(
            f"Torsion {i+1}",
            fontsize=12,
            fontweight="bold",
            pad=12
        )

        ax.grid(alpha=0.6)

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.tight_layout(rect=[0,0,1,0.95])

    _save(fig, outdir, "L_Torsions")
# ======================================================================
# 6. PROTEIN-LIGAND CONTACTS (Stacked Interaction Fractions)
# ======================================================================

def plot_contacts(indir, outdir, nframes, min_fraction=0.05):

    interaction_types = [
        ("H-bond",       "PL-Contacts_HBond.dat",       "#2e8b57"),
        ("Hydrophobic",  "PL-Contacts_Hydrophobic.dat", "#8856a7"),
        ("Ionic",        "PL-Contacts_Ionic.dat",       "#e7298a"),
        ("Water Bridge", "PL-Contacts_WaterBridge.dat", "#4eb3d3"),
        ("Pi-Cation",    "PL-Contacts_Pi-Cation.dat",   "#e6a010"),
        ("Pi-Pi",        "PL-Contacts_Pi-Pi.dat",       "#d7301f")
    ]

    fractions = {}
    residue_names = {}
    active_types = []

    for label, filename, color in interaction_types:

        filepath = os.path.join(indir, filename)

        if not os.path.exists(filepath):
            continue

        counts = Counter()

        with open(filepath) as fh:

            next(fh)

            for line in fh:

                cols = line.split()

                if len(cols) < 4:
                    continue

                residue = int(cols[1])

                counts[residue] += 1

                residue_names[residue] = cols[3]

        if not counts:
            continue

        active_types.append((label, color))

        for residue, count in counts.items():

            fractions.setdefault(residue, {})

            fractions[residue][label] = count / nframes

    if not fractions:
        print("  [skip] no contact data found")
        return

    total_fraction = {
        residue: sum(values.values())
        for residue, values in fractions.items()
    }

    residues = sorted([
        residue
        for residue, total in total_fraction.items()
        if total >= min_fraction
    ])

    labels = [
        f"{residue_names[r]}{r}"
        for r in residues
    ]

    fig_width = max(9, 0.45 * len(residues))

    fig, ax = plt.subplots(figsize=(fig_width, 5.5))

    bottom = np.zeros(len(residues))

    for label, color in active_types:

        values = np.array([
            fractions[r].get(label, 0)
            for r in residues
        ])

        ax.bar(
            range(len(residues)),
            values,
            bottom=bottom,
            color=color,
            width=0.75,
            edgecolor="white",
            linewidth=0.4,
            label=label
        )

        bottom += values

    ax.set_title(
        "Protein–Ligand Interaction Fractions",
        fontsize=18,
        fontweight="bold",
        pad=16
    )

    ax.set_ylabel(
        "Interaction Fraction",
        fontsize=15
    )

    ax.set_xlabel(
        "Protein Residue",
        fontsize=15
    )

    ax.set_xticks(range(len(residues)))

    ax.set_xticklabels(
        labels,
        rotation=90,
        fontsize=10
    )

    ax.set_xlim(-0.6, len(residues)-0.4)

    ax.set_ylim(0, bottom.max()*1.12)

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=12,
        length=6,
        width=1.2
    )

    ax.tick_params(
        axis="both",
        which="minor",
        length=3,
        width=1.0
    )

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.legend(
        loc="upper left",
        fontsize=11,
        frameon=False,
        ncol=2
    )

    fig.tight_layout()

    _save(fig, outdir, "PL_Contacts")

# ======================================================================
# main
# ======================================================================
def main():
    ap = argparse.ArgumentParser(description="Generate all SID publication figures.")
    ap.add_argument("--indir", default=".", help="folder with the .dat files")
    ap.add_argument("--outdir", default=None, help="output folder (default: --indir)")
    ap.add_argument("--time", type=float, default=200.0,
                    help="total simulation length in ns (default 200)")
    ap.add_argument("--min-fraction", type=float, default=0.05,
                    help="min total interaction fraction to show a residue")
    args = ap.parse_args()

    indir = args.indir
    outdir = args.outdir or indir
    os.makedirs(outdir, exist_ok=True)

    nframes = detect_nframes(indir)
    if nframes is None:
        nframes = int(round(args.time / 0.2))  # fallback: assume 0.2 ns/frame
        print(f"Frame count not detected; assuming {nframes}.")
    print(f"Input : {indir}\nOutput: {outdir}\nFrames: {nframes}  |  Time: {args.time} ns\n")

    plot_rmsd(indir, outdir, args.time)
    plot_protein_rmsf(indir, outdir)
    plot_ligand_rmsf(indir, outdir)
    plot_ligand_properties(indir, outdir, args.time)
    plot_ligand_torsions(indir, outdir)
    plot_contacts(indir, outdir, nframes, args.min_fraction)
    print("\nDone.")


if __name__ == "__main__":
    main()
