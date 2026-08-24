#!/usr/bin/env python3
"""
Publication-quality Protein RMSD plot (Apo MD Simulation)

Input : PL_RMSD.dat
Columns:
Frame  Prot_CA  Prot_Backbone  Prot_Sidechain  Prot_All_Heavy

Output:
P_RMSD.png
P_RMSD.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# ----------------------------------------------------
# USER SETTINGS
# ----------------------------------------------------
INFILE = "P_RMSD.dat"
OUTNAME = "P_RMSD"

TOTAL_TIME = 200.0      # ns
UNITS = "Angstrom"      # or "nm"

# Choose which RMSD to plot
PROTEIN_COL = "Prot_CA"
# Options:
# Prot_CA
# Prot_Backbone
# Prot_Sidechain
# Prot_All_Heavy

# ----------------------------------------------------
# LOAD DATA
# ----------------------------------------------------
COLNAMES = [
    "frame",
    "Prot_CA",
    "Prot_Backbone",
    "Prot_Sidechain",
    "Prot_All_Heavy"
]

col = {name:i for i,name in enumerate(COLNAMES)}

data = np.loadtxt(INFILE, comments="#")

frame = data[:, col["frame"]]
time = frame/frame[-1]*TOTAL_TIME

rmsd = data[:, col[PROTEIN_COL]]

unit = "Å" if UNITS=="Angstrom" else "nm"

# ----------------------------------------------------
# STYLE
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

COLOR = "#1f4e79"

fig, ax = plt.subplots(figsize=(7.2,4.4))

ax.plot(
    time,
    rmsd,
    color=COLOR,
    lw=1.2,
    label=PROTEIN_COL.replace("Prot_","")
)

ax.set_xlabel("Time (ns)")
ax.set_ylabel(f"Protein RMSD ({unit})")

ax.set_xlim(0,TOTAL_TIME)
ax.set_ylim(0,rmsd.max()*1.15)

ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_minor_locator(AutoMinorLocator(2))

ax.legend(frameon=False,loc="upper left")

# Mean ± SD after equilibration
eq = time >= 0.25*TOTAL_TIME

text = (
    f"Mean = {rmsd[eq].mean():.2f} ± {rmsd[eq].std():.2f} {unit}"
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

fig.savefig(f"{OUTNAME}.png",bbox_inches="tight")
fig.savefig(f"{OUTNAME}.pdf",bbox_inches="tight")

print(f"Saved {OUTNAME}.png and {OUTNAME}.pdf")
