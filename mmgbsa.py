import pandas as pd

# ==========================================================
# Read MM/GBSA results
# ==========================================================
df = pd.read_csv("trajectory-mmgbsa-prime-out.csv")

# ==========================================================
# Energy terms to summarize
# ==========================================================
energy_cols = {
    "ΔGbind": "r_psp_MMGBSA_dG_Bind",
    "ΔGbindLipo": "r_psp_MMGBSA_dG_Bind_Lipo",
    "ΔGbindvdW": "r_psp_MMGBSA_dG_Bind_vdW",
    "ΔGbindCoulomb": "r_psp_MMGBSA_dG_Bind_Coulomb",
    "ΔGbindHbond": "r_psp_MMGBSA_dG_Bind_Hbond",
    "ΔGbindSolvGB": "r_psp_MMGBSA_dG_Bind_Solv_GB",
    "ΔGbindCovalent": "r_psp_MMGBSA_dG_Bind_Covalent"
}

# ==========================================================
# Calculate Mean ± SD
# ==========================================================
summary = []

for energy_name, column in energy_cols.items():
    summary.append({
        "Energy (kcal/mol)": energy_name,
        "Mean": round(df[column].mean(), 2),
        "SD": round(df[column].std(), 2),
        "Result": f"{df[column].mean():.2f} ± {df[column].std():.2f}"
    })

summary_df = pd.DataFrame(summary)

# ==========================================================
# Save to Excel
# ==========================================================
output_file = "MMGBSA_Summary.xlsx"

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    summary_df.to_excel(writer, index=False, sheet_name="MMGBSA Summary")

print(f"Summary saved as: {output_file}")

# Display table
print(summary_df)