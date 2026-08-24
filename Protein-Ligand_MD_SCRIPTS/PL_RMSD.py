#!/usr/bin/env python3
"""
Publication-quality Protein-Ligand RMSD plot from an MD simulation.

Input : PL_RMSD.dat  (columns: frame, Prot_CA, Prot_Backbone, Prot_Sidechain,
                       Prot_All_Heavy, Lig_wrt_Protein, Lig_wrt_Ligand)
Output: PL_RMSD.png / PL_RMSD.pdf

Notes
-----
* X-axis is converted from frame index to simulation time in ns.
* RMSD values are assumed to be in Angstrom (Desmond/Maestro SID export).
  If your values are in nm (typical raw GROMACS gmx rms output), set
  UNITS = "nm" below.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# ----------------------------------------------------------------------
# 1. USER SETTINGS  -- edit these
# ----------------------------------------------------------------------
INFILE      = "PL_RMSD.dat"   # path to your data file
OUTNAME     = "PL_RMSD"       # output filename (no extension)
TOTAL_TIME  = 200.0           # total simulation length in ns
UNITS       = "Angstrom"      # "Angstrom" or "nm"

# Which protein trace to plot on the LEFT axis (pick the column name):
#   Prot_CA | Prot_Backbone | Prot_Sidechain | Prot_All_Heavy
PROTEIN_COL = "Prot_CA"
# Which ligand trace on the RIGHT axis:
#   Lig_wrt_Protein  (ligand fit on protein - stability in the pocket)
#   Lig_wrt_Ligand   (ligand internal fluctuation)
LIGAND_COL  = "Lig_wrt_Protein"

# ----------------------------------------------------------------------
# 2. LOAD DATA
# ----------------------------------------------------------------------
COLNAMES = ["frame", "Prot_CA", "Prot_Backbone", "Prot_Sidechain",
            "Prot_All_Heavy", "Lig_wrt_Protein", "Lig_wrt_Ligand"]
col = {name: i for i, name in enumerate(COLNAMES)}

data  = np.loadtxt(INFILE, comments="#")
frame = data[:, col["frame"]]

# Convert frame -> time (ns). Frame 0 = 0 ns, last frame = TOTAL_TIME ns.
time = frame / frame[-1] * TOTAL_TIME

prot = data[:, col[PROTEIN_COL]]
lig  = data[:, col[LIGAND_COL]]

unit_label = "\u00c5" if UNITS == "Angstrom" else "nm"

# ----------------------------------------------------------------------
# 3. STYLE
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "DejaVu Sans"],
    "font.size":        12,
    "axes.linewidth":   1.1,
    "xtick.direction":  "in",
    "ytick.direction":  "in",
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "savefig.dpi":      300,
    "figure.dpi":       120,
})

C_PROT = "#1f4e79"   # deep blue  - protein
C_LIG  = "#c0392b"   # brick red  - ligand

fig, ax1 = plt.subplots(figsize=(7.2, 4.4))

# ---- Protein (left axis) ----
ax1.plot(time, prot, color=C_PROT, lw=1.0, alpha=0.9,
         label=f"Protein ({PROTEIN_COL.replace('Prot_', 'C')
                          if PROTEIN_COL=='Prot_CA' else PROTEIN_COL.replace('Prot_','')})")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel(f"Protein RMSD ({unit_label})", color=C_PROT)
ax1.tick_params(axis="y", colors=C_PROT)
ax1.spines["left"].set_color(C_PROT)

# ---- Ligand (right axis) ----
ax2 = ax1.twinx()
ax2.plot(time, lig, color=C_LIG, lw=1.0, alpha=0.9,
         label=f"Ligand ({LIGAND_COL.replace('Lig_wrt_', 'wrt ')})")
ax2.set_ylabel(f"Ligand RMSD ({unit_label})", color=C_LIG)
ax2.tick_params(axis="y", colors=C_LIG)
ax2.spines["right"].set_color(C_LIG)
ax2.spines["left"].set_visible(False)

# ---- Axis limits & ticks ----
ax1.set_xlim(0, TOTAL_TIME)
ax1.xaxis.set_major_locator(MultipleLocator(50))
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))

ymax = max(prot.max(), lig.max()) * 1.15
ax1.set_ylim(0, ymax)
ax2.set_ylim(0, ymax)          # shared scale keeps the comparison honest
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))

# ---- Legend (combined) ----
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False,
           fontsize=10, handlelength=1.6)

# ---- Mean +/- SD annotation (equilibrated portion: last 75%) ----
eq = time >= 0.25 * TOTAL_TIME
txt = (f"Protein: {prot[eq].mean():.2f} \u00b1 {prot[eq].std():.2f} {unit_label}\n"
       f"Ligand:  {lig[eq].mean():.2f} \u00b1 {lig[eq].std():.2f} {unit_label}")
ax1.text(0.985, 0.04, txt, transform=ax1.transAxes, ha="right", va="bottom",
         fontsize=8.5, color="#333333",
         bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9))

fig.tight_layout()
fig.savefig(f"{OUTNAME}.png", bbox_inches="tight")
fig.savefig(f"{OUTNAME}.pdf", bbox_inches="tight")   # vector version for journals
print(f"Saved {OUTNAME}.png and {OUTNAME}.pdf")
