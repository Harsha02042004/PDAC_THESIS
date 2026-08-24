import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# MM/GBSA ENERGY COMPONENTS
# AKT1, MMP13 and P4HA3 = last 50 frames
# MEK1 = last 200 frames
# ==========================================================

components = [
    "Lipophilic",
    "van der Waals",
    "Coulomb",
    "H-bond",
    "Solvation",
    "Covalent"
]

# ==========================================================
# MEAN ENERGY VALUES
# ==========================================================

AKT1 = [
    -21.83,
    -54.90,
    -30.86,
    -2.76,
    37.17,
    3.52
]

MEK1 = [
    -15.96,
    -41.10,
    -273.55,
    -2.97,
    326.49,
    5.28
]

MMP13 = [
    -29.78,
    -85.42,
    -76.60,
    -3.30,
    99.57,
    6.26
]

P4HA3 = [
    -22.81,
    -52.10,
    -38.55,
    -3.88,
    42.58,
    2.83
]

# ==========================================================
# STANDARD DEVIATIONS
# ==========================================================

AKT1_sd = [
    1.05,
    1.84,
    4.33,
    0.45,
    2.84,
    1.80
]

MEK1_sd = [
    1.20,
    3.01,
    12.20,
    0.47,
    11.90,
    1.68
]

MMP13_sd = [
    2.18,
    7.64,
    15.62,
    0.64,
    16.05,
    6.94
]

P4HA3_sd = [
    1.32,
    1.72,
    6.53,
    0.50,
    3.65,
    1.92
]

# ==========================================================
# PLOT FIGURE 9
# ==========================================================

x = np.arange(len(components))
width = 0.20

fig, ax = plt.subplots(figsize=(12, 7))

ax.bar(
    x - 1.5 * width,
    AKT1,
    width,
    yerr=AKT1_sd,
    capsize=4,
    label="AKT1–Silibinin"
)

ax.bar(
    x - 0.5 * width,
    MEK1,
    width,
    yerr=MEK1_sd,
    capsize=4,
    label="MEK1–Silibinin"
)

ax.bar(
    x + 0.5 * width,
    MMP13,
    width,
    yerr=MMP13_sd,
    capsize=4,
    label="MMP13–Silibinin"
)

ax.bar(
    x + 1.5 * width,
    P4HA3,
    width,
    yerr=P4HA3_sd,
    capsize=4,
    label="P4HA3–Silibinin"
)

# Zero reference line
ax.axhline(
    0,
    color="black",
    linestyle="--",
    linewidth=1
)

# ==========================================================
# AXIS LABELS
# ==========================================================

ax.set_xticks(x)
ax.set_xticklabels(
    components,
    fontsize=12
)

ax.tick_params(
    axis="y",
    labelsize=12
)

ax.set_ylabel(
    "Energy (kcal/mol)",
    fontsize=14
)

ax.set_xlabel(
    "MM/GBSA Energy Components",
    fontsize=14
)

ax.set_title(
    "Comparison of MM/GBSA Energy Components",
    fontsize=16
)

# Legend
ax.legend(frameon=True)

plt.tight_layout()

# ==========================================================
# SAVE FIGURE
# ==========================================================

plt.savefig(
    "Figure9_MMGBSA_Energy_Components.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "Figure9_MMGBSA_Energy_Components.pdf",
    bbox_inches="tight"
)

plt.show()