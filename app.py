
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
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

st.dataframe(
    cln6_data.head(),
    use_container_width=True,
    hide_index=True
)

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
# Label biologically relevant genes (if significant)
important_genes = [
    "CLN6",
    "TFEB",
    "SQSTM1",
    "ATG5",
    "LAMP1",
    "CTSD"
]

for _, row in volcano_data.iterrows():
    if (
        row["Gene Symbol"] in important_genes
        and row["Significant"]
    ):
        ax.text(
            row["Difference"],
            row["minus_log10_p"],
            row["Gene Symbol"],
            fontsize=8
        )

# خطوط مرجعية
ax.axvline(-0.5, linestyle="--", color="steelblue")
ax.axvline(0.5, linestyle="--", color="steelblue")
ax.axhline(-np.log10(0.05), linestyle="--", color="darkgreen")

ax.set_xlabel("Difference")
ax.set_ylabel("-log10(P-value)")
ax.set_title("Volcano Plot of Differential Gene Expression")

ax.legend()
plt.tight_layout()
st.pyplot(fig)
# -----------------------------
# Horizontal Divider
# -----------------------------
st.divider()

# -----------------------------
# Heatmap
# -----------------------------
st.header("🔥 Heatmap of Top Differentially Expressed Genes")
# اختيار أهم 30 جين بناءً على قيمة P-value
top_genes = (
    cln6_data
    .sort_values("P_value")
    .head(30)
)

# اختيار أعمدة العينات
heatmap_data = top_genes[
    ["WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3"]
]

# جعل أسماء الجينات هي الفهرس
heatmap_data.index = top_genes["Gene Symbol"]

fig, ax = plt.subplots(figsize=(8,8))

im = ax.imshow(
    heatmap_data,
    aspect="auto",
    cmap="coolwarm"
)

# أسماء العينات
ax.set_xticks(range(len(heatmap_data.columns)))
ax.set_xticklabels(heatmap_data.columns, rotation=45)

# أسماء الجينات
ax.set_yticks(range(len(heatmap_data.index)))
ax.set_yticklabels(heatmap_data.index, fontsize=8)

# عنوان
ax.set_title(
    "Top 30 Differentially Expressed Genes",
    fontsize=14,
    fontweight="bold"
)

# Colorbar
cbar = plt.colorbar(im)
cbar.set_label("Expression Level")

plt.tight_layout()

st.pyplot(fig)

# -----------------------------
# Genes of Biological Interest
# -----------------------------
st.header("🧬 Genes of Biological Interest")

st.markdown("""
The following genes were selected based on their biological relevance
to CLN6 disease, lysosomal function, and autophagy pathways.
""")

gene_info = pd.DataFrame({
    "Gene": [
        "CLN6",
        "TFEB",
        "SQSTM1",
        "LAMP1",
        "CTSD",
        "ATG5"
    ],
    "Biological Role": [
        "Disease-causing gene",
        "Master regulator of lysosomal biogenesis",
        "Autophagy receptor",
        "Lysosomal membrane protein",
        "Lysosomal protease",
        "Autophagy-related protein"
    ]
})

st.dataframe(
    gene_info,
    use_container_width=True,
    hide_index=True
)

st.divider()
# -----------------------------
# Gene Search
# -----------------------------
st.header("🔍 Gene Search")

gene_name = st.text_input(
    "Or type a Gene Symbol",
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


st.markdown("""
---
CLN6 Epigenomic Analysis Dashboard

Developed by Sara Saad

Bachelor of Chemistry | Bioinformatics Enthusiast

Data Source: NCBI Gene Expression Omnibus (GEO)

Built using Python • Streamlit • Pandas • Matplotlib
""")
