import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ast
from pathlib import Path


# ========== CONFIGURATION ==========
DATA_PATH = Path("data/eigenfold_summary.csv")
OUTPUT_DIR = Path("data/parsed")
OUTPUT_DIR.mkdir(exist_ok=True)
sns.set(style="whitegrid", context="paper", palette="muted")


# ========== LOAD & CLEAN ==========
print("[INFO] Loading data from:", DATA_PATH)
df = pd.read_csv(DATA_PATH)

# Convert stringified lists to Python lists, if needed
if df["stable_states"].dtype == object:
    try:
        df["stable_states"] = df["stable_states"].apply(ast.literal_eval)
    except Exception:
        print("[WARN] Could not parse 'stable_states' — skipping literal_eval.")

# Derived columns
df["R"] = df["R"].astype(int)
df["n_stable"] = df["n_stable"].astype(int)
df["entropy"] = df["entropy"].astype(float)
df["time"] = df["time"].astype(float)
df["dominant_p"] = df["dominant_p"].astype(float)
df["stable_fraction"] = df["stable_p_sum"].astype(float)


# ========== SUMMARY TABLES ==========
print("\n=== Global Summary ===")
summary = df.describe()[["dominant_p", "entropy", "n_stable", "energy", "time"]]
print(summary)

print("\n=== Lattice Summary ===")
lattice_summary = (
    df.groupby("lattice")[["dominant_p", "entropy", "n_stable", "energy", "time"]]
    .mean()
    .round(3)
)
print(lattice_summary)
lattice_summary.to_csv(OUTPUT_DIR / "lattice_summary.csv")


print("\n=== Scaling Summary ===")
scaling_summary = (
    df.groupby("R")[["n_qubits", "time", "entropy", "n_stable"]]
    .mean()
    .round(3)
)
print(scaling_summary)
scaling_summary.to_csv(OUTPUT_DIR / "scaling_summary.csv")


# ========== PLOTS ==========

def saveplot(fig, name: str):
    path = OUTPUT_DIR / f"{name}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Saved {path}")


# --- 1. Runtime vs R and Lattice ---
fig, ax = plt.subplots(figsize=(7, 5))
sns.lineplot(data=df, x="R", y="time", hue="lattice", style="lattice",
             markers=True, err_style="bars", ax=ax)
ax.set_title("Runtime Scaling with Peptide Length across Lattice Types")
ax.set_xlabel("Residue Count (R)")
ax.set_ylabel("Runtime (s)")
saveplot(fig, "runtime_scaling")


# --- 2. Entropy vs Lattice ---
fig, ax = plt.subplots(figsize=(6, 4))
sns.boxplot(data=df, x="lattice", y="entropy", ax=ax)
ax.set_title("Entropy Distribution across Lattice Types")
ax.set_xlabel("Lattice Type")
ax.set_ylabel("Entropy")
saveplot(fig, "entropy_by_lattice")


# --- 3. Stable States vs Lattice ---
fig, ax = plt.subplots(figsize=(6, 4))
sns.boxplot(data=df, x="lattice", y="n_stable", ax=ax)
ax.set_title("Number of Stable Conformations across Lattices")
ax.set_xlabel("Lattice")
ax.set_ylabel("Stable Conformations (p > 0.1)")
saveplot(fig, "stable_states_by_lattice")


# --- 4. Energy vs Entropy ---
fig, ax = plt.subplots(figsize=(6, 4))
sns.scatterplot(data=df, x="entropy", y="energy", hue="R", palette="viridis", ax=ax)
ax.set_title("Energy Landscape: Entropy vs Energy")
ax.set_xlabel("Entropy")
ax.set_ylabel("Energy (a.u.)")
saveplot(fig, "energy_vs_entropy")


# --- 5. Correlation Matrix ---
corr = df[["dominant_p", "entropy", "n_stable", "energy", "time"]].corr()
fig, ax = plt.subplots(figsize=(6, 5))
sns.heatmap(corr, annot=True, cmap="coolwarm", fmt=".2f", ax=ax)
ax.set_title("Correlation Matrix of EigenFold Metrics")
saveplot(fig, "correlation_matrix")


# --- 6. Lattice-wise Mean Runtime with Error Bars ---
time_summary = (
    df.groupby(["R", "lattice"])["time"]
    .agg(["mean", "std"])
    .reset_index()
)
fig, ax = plt.subplots(figsize=(7, 5))
for lattice in df["lattice"].unique():
    sub = time_summary[time_summary["lattice"] == lattice]
    ax.errorbar(sub["R"], sub["mean"], yerr=sub["std"], marker="o", label=lattice.upper())
ax.set_xlabel("Residue Count (R)")
ax.set_ylabel("Runtime (s)")
ax.set_title("Runtime Scaling (mean ± std)")
ax.legend()
saveplot(fig, "runtime_errorbars")


# --- 7. Fraction of Multi-Stable Systems ---
multi_frac = (df["n_stable"] > 1).mean() * 100
print(f"\n[INFO] {multi_frac:.1f}% of systems exhibit multiple stable conformations (p > 0.1).")

fig, ax = plt.subplots(figsize=(6, 4))
sns.barplot(data=df, x="lattice", y="stable_fraction", ax=ax)
ax.set_title("Fraction of Probability Mass in Stable Conformations")
ax.set_xlabel("Lattice Type")
ax.set_ylabel("Stable Probability Fraction")
saveplot(fig, "stable_fraction_by_lattice")


# --- 8. Optional: Runtime vs Energy Scatter ---
fig, ax = plt.subplots(figsize=(6, 4))
sns.scatterplot(data=df, x="energy", y="time", hue="lattice", style="R", ax=ax)
ax.set_title("Runtime vs Minimum Energy by Lattice Type")
ax.set_xlabel("Minimum Energy (a.u.)")
ax.set_ylabel("Runtime (s)")
saveplot(fig, "runtime_vs_energy")


# --- 9. Summary report ---
report = {
    "num_systems": len(df),
    "mean_runtime": df["time"].mean(),
    "mean_entropy": df["entropy"].mean(),
    "mean_stable_states": df["n_stable"].mean(),
    "multi_stable_fraction": multi_frac,
}
print("\n=== Summary Report ===")
for k, v in report.items():
    print(f"{k:>20}: {v:.3f}")

(pd.Series(report).to_frame("value")
 .to_csv(OUTPUT_DIR / "summary_report.csv"))

print("\n[INFO] Analysis complete. All figures saved in:", OUTPUT_DIR)