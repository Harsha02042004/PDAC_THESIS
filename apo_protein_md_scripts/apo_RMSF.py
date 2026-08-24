#!/usr/bin/env python3
"""
Publication-quality Protein RMSF plot (Apo Simulation)

Input:
    P_RMSF.dat

Output:
    P_RMSF.png
    P_RMSF.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# ----------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------
INFILE = "P_RMSF.dat"
OUTNAME = "P_RMSF"

# Choose which RMSF column to plot
RMSF_COL = "CA"

# Options:
# CA
# Backbone
# Sidechain
# All_Heavy
# B-factor

UNIT = "Å"

# ----------------------------------------------------
# Column indices in Desmond SID export
# ----------------------------------------------------
COLIDX = {
    "CA": 4,
    "Backbone": 5,
    "Sidechain": 6,
    "All_Heavy": 7,
    "B-factor": 8
}

# ----------------------------------------------------
# Read data
# ----------------------------------------------------
residue = []
values = []

with open(INFILE) as f:
    next(f)          # Skip header

    for line in f:

        parts = line.split()

        if len(parts) < 9:
            continue

        residue.append(int(parts[0]))
        values.append(float(parts[COLIDX[RMSF_COL]]))

residue = np.array(residue)
values = np.array(values)

# ----------------------------------------------------
# Plot style
# ----------------------------------------------------
plt.rcParams.update({

    "font.family":"sans-serif",
    "font.sans-serif":["Arial","DejaVu Sans"],

    "font.size":12,
    "axes.linewidth":1.2,

    "xtick.direction":"in",
    "ytick.direction":"in",

    "xtick.major.size":5,
    "ytick.major.size":5,

    "xtick.minor.size":3,
    "ytick.minor.size":3,

    "savefig.dpi":300,
    "figure.dpi":120
})

COLOR = "#5b7488"

fig, ax = plt.subplots(figsize=(8,4.5))

ax.plot(
    residue,
    values,
    color=COLOR,
    linewidth=1.8
)

ax.set_title("Protein RMSF")

ax.set_xlabel("Residue Index")
ax.set_ylabel(f"RMSF ({UNIT})")

ax.set_xlim(residue.min(), residue.max())
ax.set_ylim(0, values.max()*1.10)

ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_minor_locator(AutoMinorLocator(2))

ax.legend(
    [RMSF_COL if RMSF_COL != "CA" else "Cα"],
    frameon=False,
    loc="upper right"
)

# Mean ± SD annotation
text = (
    f"Mean = {values.mean():.2f} ± {values.std():.2f} {UNIT}"
)

ax.text(
    0.985,
    0.04,
    text,
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=9,
    bbox=dict(
        boxstyle="round,pad=0.4",
        fc="white",
        ec="#cccccc",
        alpha=0.9
    )
)

fig.tight_layout()

fig.savefig(f"{OUTNAME}.png", bbox_inches="tight")
fig.savefig(f"{OUTNAME}.pdf", bbox_inches="tight")

print(f"Saved {OUTNAME}.png and {OUTNAME}.pdf")
