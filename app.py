
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Page Configuration
# -----------------------------
st.set_page_config(
    page_title="CLN6 Epigenomic Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# -----------------------------
# Dashboard Title
# -----------------------------
st.title("🧬 CLN6 Epigenomic Analysis Dashboard")

st.markdown("""
This interactive dashboard explores differential gene expression
associated with CLN6 mutation using a publicly available processed
microarray dataset from the NCBI Gene Expression Omnibus (GEO).

The dashboard provides:

- 📊 Dataset overview
- 🌋 Volcano Plot for differential expression
- 🔍 Gene-level exploration
- 📈 Summary statistics

Purpose: Demonstrate bioinformatics data exploration and visualization
using Python and Streamlit.
""")

# -----------------------------
# Load Dataset
# -----------------------------
try:
    cln6_data = pd.read_csv("CLN6_dashboard_data.csv")
except FileNotFoundError:
    st.error("Dataset not found. Please upload CLN6_dashboard_data.csv")
    st.stop()
    
# -----------------------------
# About CLN6
# -----------------------------
st.subheader("🧬 About CLN6")

st.markdown("""
CLN6 is a rare neurodegenerative disorder belonging to the
Neuronal Ceroid Lipofuscinoses (NCLs).

Mutations in the CLN6 gene lead to progressive neuronal dysfunction
through complex molecular mechanisms.

As research continues, understanding this disease requires connecting
multiple biological pathways rather than focusing on a single mechanism.

This dashboard provides an interactive way to explore gene expression
changes and compare wild-type and CLN6 mutant samples using a
publicly available processed microarray dataset.
""")

st.divider()

# -----------------------------
# Project Summary
# -----------------------------
st.header("📌 Project Summary")

col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        label="🧪 Samples",
        value="6"
    )

with col2:
    st.metric(
        label="👥 Experimental Groups",
        value="3 WT | 3 Mut"
    )

with col3:
    st.metric(
        label="🧬 Genes Analysed",
        value=f"{len(cln6_data):,}"
    )
  
st.divider()

# -----------------------------
# Dataset Information
# -----------------------------

st.subheader("📂 Dataset Information")

col1, col2 = st.columns(2)

with col1:
    st.info("""
Source

NCBI Gene Expression Omnibus (GEO)

Processed Microarray Dataset
""")

with col2:
    st.info("""
Study Design

3 Wild Type (WT)

3 CLN6 Mutant
""")

st.divider()

# -----------------------------
# Sidebar
# -----------------------------
st.sidebar.title("🧬 Project Information")

st.sidebar.markdown("""
### CLN6 Epigenomic Analysis

Disease: Neuronal Ceroid Lipofuscinosis Type 6 (CLN6)

Dataset Source: NCBI GEO

Platform: Microarray Gene Expression

Samples: 3 Wild Type (WT) + 3 Mutant (Mut)

Objective:
Explore gene expression changes associated with CLN6 mutation and identify genes potentially involved in early disease mechanisms.
""")

# -----------------------------
# Dataset Preview
# -----------------------------
st.header("📋 Dataset Preview")

st.markdown(
    "Browse the processed CLN6 gene expression dataset used in this project."
)

st.write(cln6_data.head())

st.divider()
# -----------------------------
# Quick Statistics
# -----------------------------
st.header("📊 Dataset Statistics")

col1, col2, col3 = st.columns(3)

with col1:
    st.metric("Total Genes", f"{len(cln6_data):,}")

with col2:
    significant = (cln6_data["P_value"] < 0.05).sum()
    st.metric("Significant Genes", f"{significant:,}")

with col3:
    avg_difference = cln6_data["Difference"].abs().mean()
    st.metric("Average |Difference|", f"{avg_difference:.3f}")
st.divider()
# -----------------------------
# Volcano Plot
# -----------------------------
import numpy as np

st.header("🌋 Volcano Plot")

volcano_data = cln6_data.copy()

# حساب -log10(P-value)
volcano_data["minus_log10_p"] = -np.log10(
    volcano_data["P_value"].replace(0, 1e-300)
)

# تحديد الجينات المهمة
volcano_data["Significant"] = (
    (volcano_data["P_value"] < 0.05) &
    (volcano_data["Difference"].abs() > 0.5)
)

fig, ax = plt.subplots(figsize=(8,6))

# الجينات غير المهمة
ax.scatter(
    volcano_data.loc[~volcano_data["Significant"], "Difference"],
    volcano_data.loc[~volcano_data["Significant"], "minus_log10_p"],
    color="lightgray",
    alpha=0.6,
    s=20,
    label="Not Significant"
)

# الجينات المهمة
ax.scatter(
    volcano_data.loc[volcano_data["Significant"], "Difference"],
    volcano_data.loc[volcano_data["Significant"], "minus_log10_p"],
    color="crimson",
    alpha=0.8,
    s=25,
    label="Significant"
)

# خطوط مرجعية
ax.axvline(-0.5, linestyle="--", color="steelblue")
ax.axvline(0.5, linestyle="--", color="steelblue")
ax.axhline(-np.log10(0.05), linestyle="--", color="darkgreen")

ax.set_xlabel("Difference")
ax.set_ylabel("-log10(P-value)")
ax.set_title("Volcano Plot of Differential Gene Expression")

ax.legend()

st.pyplot(fig)
# -----------------------------
# Horizontal Divider
# -----------------------------
st.divider()

# -----------------------------
# Data Exploration
# -----------------------------
st.header("🔍 Data Exploration")

st.markdown("""
Explore the processed gene expression dataset and visualize
the differential expression results.
""")
# -----------------------------
# Prepare Volcano Plot Data
# -----------------------------
import numpy as np

volcano_data = cln6_data.copy()

# Calculate -log10(P-value)
volcano_data["minus_log10_p"] = -np.log10(
    volcano_data["P_value"].clip(lower=1e-300)
)

# Define significant genes
volcano_data["Significant"] = (
    (volcano_data["P_value"] < 0.05) &
    (volcano_data["Difference"].abs() > 0.5)
)

 # -----------------------------
# Volcano Plot Figure
# -----------------------------
fig, ax = plt.subplots(figsize=(9,7))

# Non-significant genes
ax.scatter(
    volcano_data.loc[~volcano_data["Significant"], "Difference"],
    volcano_data.loc[~volcano_data["Significant"], "minus_log10_p"],
    color="lightgray",
    alpha=0.6,
    s=18,
    label="Not Significant"
)

# Significant genes
ax.scatter(
    volcano_data.loc[volcano_data["Significant"], "Difference"],
    volcano_data.loc[volcano_data["Significant"], "minus_log10_p"],
    color="crimson",
    alpha=0.85,
    s=22,
    label="Significant"
)
# -----------------------------
# Volcano Plot Styling
# -----------------------------

# Threshold lines
ax.axvline(-0.5, color="steelblue", linestyle="--", linewidth=1)
ax.axvline(0.5, color="steelblue", linestyle="--", linewidth=1)
ax.axhline(-np.log10(0.05), color="darkgreen", linestyle="--", linewidth=1)

# Labels
ax.set_title(
    "Volcano Plot of Differential Gene Expression",
    fontsize=14,
    fontweight="bold"
)

ax.set_xlabel("Mean Expression Difference (Mut - WT)")
ax.set_ylabel("-log10(P-value)")

# Clean appearance
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False)

plt.tight_layout()

st.pyplot(fig)


# -----------------------------
# Gene Search
# -----------------------------
st.header("🧬 Gene Search")

gene_name = st.text_input(
    "Enter Gene Symbol",
    placeholder="Example: TFEB"
)

if gene_name:

    result = cln6_data[
        cln6_data["Gene Symbol"]
        .astype(str)
        .str.upper()
        == gene_name.upper()
    ]

    if result.empty:
        st.warning("Gene not found.")
    else:
        st.success(f"Found {len(result)} matching gene(s).")
        st.dataframe(result, use_container_width=True, hide_index=True)

# -----------------------------
# About CLN6
# -----------------------------
st.subheader("🧬 About CLN6")

st.markdown("""
CLN6 is a rare neurodegenerative disorder belonging to the
Neuronal Ceroid Lipofuscinoses (NCLs).

Mutations in the CLN6 gene lead to progressive neuronal dysfunction
through complex molecular mechanisms.

As research continues, understanding this disease requires connecting
multiple biological pathways rather than focusing on a single mechanism.

This dashboard provides an interactive way to explore gene expression
changes and compare wild-type and CLN6 mutant samples using a
publicly available processed microarray dataset.
""")

st.divider()
