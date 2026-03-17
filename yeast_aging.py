#!/usr/bin/env python3
"""
Information‑theoretic analysis of yeast scRNA‑seq data.
Computes per‑cell entropy, compares high vs low entropy cells,
and analyses changes over time. Generates a comprehensive PDF report.
"""

import os
import logging
import argparse
from pathlib import Path
import pickle
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import f_oneway, kruskal, ttest_ind
import xml.etree.ElementTree as ET
from gprofiler import GProfiler

from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch

# ---------------------------- CONFIGURATION ----------------------------
# Paths (adjust as needed)
DATA_DIR = Path("/home/anirudh/yeast_aging")
EXPR_FILE = DATA_DIR / "data/GSE210032_Log2NormCounts_1W-A37B43C45_H2_NoERCC_CSV.csv"
META_FILE = DATA_DIR / "GSE210032_family.xml"
OUTPUT_DIR = Path("outputsGSE210032")
CACHE_DIR = OUTPUT_DIR / "cache"
TOP_N_GENES = 50          # number of top genes for GO / overlap
PVAL_THRESH = 0.05
LOG2FC_THRESH = 1.0


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# Cytoskeleton gene symbols (S. cerevisiae)
CYTO_SYMBOLS = [
    "ACT1", "ABP1", "SAC6", "COF1", "TPM1", "TPM2", "MYO1", "MYO2",
    "BNI1", "BNR1", "ARP2", "ARP3", "LAS17", "SLA1", "SLA2",
    "CDC42", "RHO1", "RHO3"
]

# ---------------------------- SETUP LOGGING ----------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ---------------------------- HELPER FUNCTIONS ----------------------------
def shannon_entropy(x):
    """Compute Shannon entropy of a count vector (linear scale)."""
    x = x[x > 0]
    if len(x) == 0:
        return 0.0
    p = x / x.sum()
    return -np.sum(p * np.log2(p))

def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)

def load_expression(file_path):
    """Load log2-normalised counts and convert to linear scale."""
    logger.info(f"Loading expression data from {file_path}")
    expr_log2 = pd.read_csv(file_path, index_col=0)
    expr_lin = np.exp2(expr_log2) - 1
    logger.info(f"Loaded {expr_lin.shape[0]} genes and {expr_lin.shape[1]} cells")
    return expr_lin

def parse_metadata(file_path):
    """Parse GEO XML metadata to extract sample time points."""
    logger.info(f"Parsing metadata from {file_path}")
    tree = ET.parse(file_path)
    root = tree.getroot()
    ns = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
    records = []
    for sample in root.findall(".//geo:Sample", ns):
        title = sample.findtext("geo:Title", namespaces=ns)
        acc = sample.findtext("geo:Accession", namespaces=ns)
        time = None
        for ch in sample.findall(".//geo:Characteristics", ns):
            if "time" in ch.attrib.get("tag", "").lower():
                time = ch.text.strip()
        records.append({"accession": acc, "title": title, "time": time})
    meta = pd.DataFrame(records)
    meta["sample"] = meta["title"].str.extract(r"sc_([A-Z]\d+)")
    meta["time_hr"] = meta["time"].str.extract(r"(\d+)").astype(float)
    meta = meta.set_index("sample")
    logger.info(f"Metadata parsed: {len(meta)} samples")
    return meta

def align_metadata_expression(meta, expr):
    """Keep only samples present in both metadata and expression."""
    common = meta.index.intersection(expr.columns)
    if len(common) < len(expr.columns):
        logger.warning(f"Dropping {len(expr.columns) - len(common)} samples not in metadata")
    meta = meta.loc[common]
    expr = expr[common]
    return meta, expr

def run_gprofiler(queries, organism="scerevisiae", **kwargs):
    """
    Run g:Profiler GO enrichment with caching.
    queries: list of gene IDs or a single list.
    Returns a DataFrame (or None if failed).
    """
    gp = GProfiler(return_dataframe=True)
    # Use a hash of the query list for cache filename
    query_hash = hash(frozenset(queries))
    cache_file = CACHE_DIR / f"gprofiler_enrich_{query_hash}.pkl"
    if cache_file.exists():
        logger.info(f"Loading cached g:Profiler enrichment from {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    try:
        logger.info("Querying g:Profiler for GO enrichment...")
        result = gp.profile(organism=organism, query=queries, **kwargs)
        with open(cache_file, "wb") as f:
            pickle.dump(result, f)
        return result
    except Exception as e:
        logger.error(f"g:Profiler enrichment failed: {e}")
        return None

def convert_genes(queries, organism="scerevisiae"):
    """
    Convert gene symbols/IDs to stable identifiers using g:Profiler convert.
    Returns a DataFrame with columns: incoming, converted, name, etc.
    Caches results.
    """
    gp = GProfiler(return_dataframe=True)
    query_hash = hash(frozenset(queries))
    cache_file = CACHE_DIR / f"gprofiler_convert_{query_hash}.pkl"
    if cache_file.exists():
        logger.info(f"Loading cached gene conversion from {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    try:
        logger.info("Querying g:Profiler for gene conversion...")
        result = gp.convert(organism=organism, query=queries)
        with open(cache_file, "wb") as f:
            pickle.dump(result, f)
        return result
    except Exception as e:
        logger.error(f"g:Profiler conversion failed: {e}")
        return None

def map_gene_names(genes):
    """Convert ORF IDs to gene symbols using g:Profiler."""
    mapping_df = convert_genes(genes)
    if mapping_df is not None and 'incoming' in mapping_df.columns and 'name' in mapping_df.columns:
        return dict(zip(mapping_df["incoming"], mapping_df["name"]))
    else:
        logger.warning("Gene mapping failed; returning empty dict")
        return {}

# ---------------------------- MAIN ANALYSIS FUNCTIONS ----------------------------
def compute_entropy(expr_lin):
    logger.info("Computing per‑cell entropy")
    return expr_lin.apply(shannon_entropy, axis=0)

def test_entropy_across_time(entropy, meta):
    times = sorted(meta["time_hr"].unique())
    groups = [entropy[meta["time_hr"] == t] for t in times]
    anova_p = f_oneway(*groups).pvalue
    kruskal_p = kruskal(*groups).pvalue
    variance = {t: np.var(entropy[meta["time_hr"] == t]) for t in times}
    logger.info(f"ANOVA p={anova_p:.3e}, Kruskal p={kruskal_p:.3e}")
    return anova_p, kruskal_p, variance, times

def differential_expression_entropy(expr_lin, entropy, top_n=10):
    """High vs low entropy cells."""
    high_cells = entropy.sort_values(ascending=False).head(top_n).index
    low_cells = entropy.sort_values().head(top_n).index
    logger.info(f"Comparing {len(high_cells)} high‑entropy vs {len(low_cells)} low‑entropy cells")

    # Vectorised t-test
    high_vals = expr_lin[high_cells]
    low_vals = expr_lin[low_cells]
    t_stat, pvals = ttest_ind(high_vals, low_vals, axis=1, equal_var=False)

    mean_high = high_vals.mean(axis=1)
    mean_low = low_vals.mean(axis=1)
    log2fc = np.log2((mean_high + 1) / (mean_low + 1))

    volcano = pd.DataFrame({
        "gene": expr_lin.index,
        "log2FC": log2fc,
        "pval": pvals,
        "-log10p": -np.log10(pvals)
    }).sort_values("pval")
    return volcano, high_cells, low_cells

def differential_expression_time(expr_lin, meta):
    """Late vs early time points."""
    times = sorted(meta["time_hr"].unique())
    early_cells = meta[meta["time_hr"] == times[0]].index
    late_cells = meta[meta["time_hr"] == times[-1]].index
    logger.info(f"Comparing {len(late_cells)} late (t={times[-1]}) vs {len(early_cells)} early (t={times[0]}) cells")

    early_vals = expr_lin[early_cells]
    late_vals = expr_lin[late_cells]
    t_stat, pvals = ttest_ind(late_vals, early_vals, axis=1, equal_var=False)

    mean_late = late_vals.mean(axis=1)
    mean_early = early_vals.mean(axis=1)
    log2fc = np.log2((mean_late + 1) / (mean_early + 1))

    volcano = pd.DataFrame({
        "gene": expr_lin.index,
        "log2FC": log2fc,
        "pval": pvals,
        "-log10p": -np.log10(pvals)
    }).sort_values("pval")
    return volcano, early_cells, late_cells

def cytoskeleton_analysis(expr_lin, meta, cyto_orfs, mapping_dict):
    """Extract cytoskeleton gene expression and compute entropy and time trends."""
    cyto_expr = expr_lin.loc[expr_lin.index.intersection(cyto_orfs)]
    cyto_entropy = cyto_expr.apply(shannon_entropy, axis=0)

    # Time trajectories
    times = sorted(meta["time_hr"].unique())
    records = []
    for gene in cyto_expr.index:
        for t in times:
            vals = cyto_expr.loc[gene, meta[meta["time_hr"] == t].index]
            records.append({
                "gene": gene,
                "gene_name": mapping_dict.get(gene, gene),
                "time": t,
                "mean_expr": vals.mean(),
                "var_expr": vals.var()
            })
    cyto_time = pd.DataFrame(records)
    return cyto_entropy, cyto_time

def classify_genes(volcano_entropy, volcano_time, all_genes, top_n=50):
    """Classify genes based on membership in top DE lists."""
    entropy_top = set(volcano_entropy.head(top_n)["gene"])
    time_top = set(volcano_time.head(top_n)["gene"])

    def classify(g):
        in_e = g in entropy_top
        in_t = g in time_top
        if in_e and in_t:
            return "core_aging"
        if in_e:
            return "entropy_specific"
        if in_t:
            return "time_specific"
        return "background"

    classes = pd.DataFrame({
        "gene": all_genes,
        "class": [classify(g) for g in all_genes]
    })
    return classes, entropy_top, time_top

# ---------------------------- PLOTTING FUNCTIONS ----------------------------
def plot_entropy_box(entropy, meta, out_path):
    plt.figure()
    sns.boxplot(x=meta["time_hr"], y=entropy.values)
    sns.stripplot(x=meta["time_hr"], y=entropy.values, color="black", alpha=0.4)
    plt.xlabel("Time (hr)")
    plt.ylabel("Entropy (bits)")
    plt.title("Global Entropy per Cell")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    logger.info(f"Saved {out_path}")

def plot_cyto_entropy(cyto_entropy, meta, out_path):
    plt.figure()
    sns.boxplot(x=meta["time_hr"], y=cyto_entropy.values)
    sns.stripplot(x=meta["time_hr"], y=cyto_entropy.values, color="black", alpha=0.4)
    plt.xlabel("Time (hr)")
    plt.ylabel("Entropy (bits)")
    plt.title("Cytoskeletal Gene Entropy")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    logger.info(f"Saved {out_path}")

def plot_variance_trend(variance_dict, times, out_path):
    plt.figure()
    plt.plot(times, list(variance_dict.values()), marker='o')
    plt.xlabel("Time (hr)")
    plt.ylabel("Variance of Entropy")
    plt.title("Entropy Variance Across Time")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    logger.info(f"Saved {out_path}")

def plot_volcano(volcano, highlight_genes=None, title="", out_path=None):
    plt.figure()
    # background
    bg = volcano[~volcano["gene"].isin(highlight_genes)] if highlight_genes is not None else volcano
    plt.scatter(bg["log2FC"], bg["-log10p"], c='blue', alpha=0.4, s=10, label='other')
    if highlight_genes is not None:
        fg = volcano[volcano["gene"].isin(highlight_genes)]
        plt.scatter(fg["log2FC"], fg["-log10p"], c='red', s=20, label='cytoskeleton')
    plt.axhline(-np.log10(PVAL_THRESH), linestyle='--', color='grey', linewidth=0.8)
    plt.axvline(LOG2FC_THRESH, linestyle='--', color='grey', linewidth=0.8)
    plt.axvline(-LOG2FC_THRESH, linestyle='--', color='grey', linewidth=0.8)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if out_path:
        plt.savefig(out_path, dpi=150)
        plt.close()
        logger.info(f"Saved {out_path}")
    else:
        plt.show()

def plot_cyto_trajectories(cyto_time, out_path):
    plt.figure(figsize=(8,6))
    for gene in cyto_time["gene"].unique():
        sub = cyto_time[cyto_time["gene"] == gene]
        plt.plot(sub["time"], sub["mean_expr"], alpha=0.4)
    plt.xlabel("Time (hr)")
    plt.ylabel("Mean Expression")
    plt.title("Cytoskeletal Gene Expression Trajectories")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    logger.info(f"Saved {out_path}")

# ---------------------------- REPORT GENERATION ----------------------------
def generate_pdf_report(output_path, results):
    """
    results: dict containing all data needed for the report.
    """
    doc = SimpleDocTemplate(str(output_path), pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    # Title
    story.append(Paragraph("Information‑Theoretic Analysis of Yeast scRNA‑seq", styles['Title']))
    story.append(Spacer(1, 12))

    # Methods
    story.append(Paragraph("Methods", styles['Heading2']))
    story.append(Paragraph(
        "Entropy was computed per cell after converting log2‑normalized data to linear scale. "
        "Cells were grouped by timepoint. Differential expression was computed for high vs low entropy "
        "cells and early vs late timepoints using Welch’s t‑test. GO enrichment was performed using g:Profiler. "
        f"Thresholds: p < {PVAL_THRESH}, |log2FC| > {LOG2FC_THRESH}.",
        styles['BodyText']
    ))
    story.append(Spacer(1, 12))

    # Global entropy statistics
    story.append(Paragraph("Global Entropy Across Time", styles['Heading2']))
    story.append(Paragraph(
        f"ANOVA p‑value: {results['anova_p']:.3e}<br/>"
        f"Kruskal–Wallis p‑value: {results['kruskal_p']:.3e}",
        styles['BodyText']
    ))
    story.append(Spacer(1, 12))

    # Entropy boxplot
    story.append(Image(str(results['entropy_box']), width=400, height=250))
    story.append(Spacer(1, 12))

    # Cytoskeleton entropy boxplot
    story.append(Paragraph("Cytoskeletal Gene Entropy", styles['Heading3']))
    story.append(Image(str(results['cyto_entropy_plot']), width=400, height=250))
    story.append(Spacer(1, 12))

    # Entropy variance
    story.append(Paragraph("Variance of Entropy Over Time", styles['Heading3']))
    story.append(Image(str(results['variance_plot']), width=400, height=250))
    story.append(Spacer(1, 12))

    # Volcano plots
    story.append(Paragraph("Differential Expression", styles['Heading2']))
    story.append(Paragraph("High‑ vs Low‑Entropy Cells", styles['Heading3']))
    story.append(Image(str(results['volcano_entropy']), width=400, height=250))
    story.append(Spacer(1, 6))

    # Top DE genes table (entropy)
    top_entropy = results['volcano_entropy_df'].head(10)[['gene_name', 'log2FC', 'pval']].copy()
    top_entropy['pval'] = top_entropy['pval'].map(lambda x: f"{x:.2e}")
    data = [['Gene', 'log2FC', 'p‑value']] + top_entropy.values.tolist()
    t = Table(data, colWidths=[1.5*inch, 1*inch, 1.5*inch])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.grey),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,0), 10),
        ('BOTTOMPADDING', (0,0), (-1,0), 8),
        ('BACKGROUND', (0,1), (-1,-1), colors.beige),
        ('GRID', (0,0), (-1,-1), 1, colors.black)
    ]))
    story.append(t)
    story.append(Spacer(1, 12))

    story.append(Paragraph("Late vs Early Time Points", styles['Heading3']))
    story.append(Image(str(results['volcano_time']), width=400, height=250))
    story.append(Spacer(1, 6))

    top_time = results['volcano_time_df'].head(10)[['gene_name', 'log2FC', 'pval']].copy()
    top_time['pval'] = top_time['pval'].map(lambda x: f"{x:.2e}")
    data = [['Gene', 'log2FC', 'p‑value']] + top_time.values.tolist()
    t = Table(data, colWidths=[1.5*inch, 1*inch, 1.5*inch])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.grey),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,0), 10),
        ('BOTTOMPADDING', (0,0), (-1,0), 8),
        ('BACKGROUND', (0,1), (-1,-1), colors.beige),
        ('GRID', (0,0), (-1,-1), 1, colors.black)
    ]))
    story.append(t)
    story.append(Spacer(1, 12))

    # GO enrichment tables
    if results['go_entropy'] is not None:
        story.append(Paragraph("GO Enrichment (High‑ vs Low‑Entropy Genes)", styles['Heading2']))
        go_top = results['go_entropy'].head(10)[['native', 'name', 'p_value']].copy()
        go_top['p_value'] = go_top['p_value'].map(lambda x: f"{x:.2e}")
        data = [['GO term', 'Description', 'p‑value']] + go_top.values.tolist()
        t = Table(data, colWidths=[1.2*inch, 2.5*inch, 1*inch])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.grey),
            ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTSIZE', (0,0), (-1,0), 9),
            ('BOTTOMPADDING', (0,0), (-1,0), 8),
            ('BACKGROUND', (0,1), (-1,-1), colors.beige),
            ('GRID', (0,0), (-1,-1), 1, colors.black)
        ]))
        story.append(t)
        story.append(Spacer(1, 12))

    if results['go_time'] is not None:
        story.append(Paragraph("GO Enrichment (Late vs Early Genes)", styles['Heading2']))
        go_top = results['go_time'].head(10)[['native', 'name', 'p_value']].copy()
        go_top['p_value'] = go_top['p_value'].map(lambda x: f"{x:.2e}")
        data = [['GO term', 'Description', 'p‑value']] + go_top.values.tolist()
        t = Table(data, colWidths=[1.2*inch, 2.5*inch, 1*inch])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.grey),
            ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTSIZE', (0,0), (-1,0), 9),
            ('BOTTOMPADDING', (0,0), (-1,0), 8),
            ('BACKGROUND', (0,1), (-1,-1), colors.beige),
            ('GRID', (0,0), (-1,-1), 1, colors.black)
        ]))
        story.append(t)
        story.append(Spacer(1, 12))

    # Cytoskeleton trajectories
    story.append(Paragraph("Cytoskeleton Gene Expression Over Time", styles['Heading2']))
    story.append(Image(str(results['cyto_trajectories']), width=400, height=250))
    story.append(Spacer(1, 12))

    # Classification summary
    story.append(Paragraph("Gene Classification", styles['Heading2']))
    class_counts = results['classification']['class'].value_counts()
    story.append(Paragraph(
        f"Core aging (both lists): {class_counts.get('core_aging', 0)}<br/>"
        f"Entropy‑specific: {class_counts.get('entropy_specific', 0)}<br/>"
        f"Time‑specific: {class_counts.get('time_specific', 0)}<br/>"
        f"Background: {class_counts.get('background', 0)}",
        styles['BodyText']
    ))
    story.append(Spacer(1, 6))

    overlap_list = results['entropy_top'].intersection(results['time_top'])
    story.append(Paragraph(f"Overlap of top {TOP_N_GENES} genes: {len(overlap_list)} genes", styles['BodyText']))
    if overlap_list:
        story.append(Paragraph("Genes in overlap: " + ", ".join(list(overlap_list)[:20]) + 
                               ("..." if len(overlap_list)>20 else ""), styles['BodyText']))

    # Build PDF
    doc.build(story)
    logger.info(f"PDF report generated: {output_path}")

# ---------------------------- MAIN ----------------------------
def main():
    parser = argparse.ArgumentParser(description="Yeast aging scRNA‑seq entropy analysis")
    parser.add_argument("--expr", type=str, default=str(EXPR_FILE), help="Expression CSV file")
    parser.add_argument("--meta", type=str, default=str(META_FILE), help="Metadata XML file")
    parser.add_argument("--out", type=str, default="outputs", help="Output directory")
    parser.add_argument("--top_n", type=int, default=TOP_N_GENES, help="Number of top genes for GO/overlap")
    args = parser.parse_args()

    # Setup directories
    out_dir = Path(args.out)
    ensure_dir(out_dir)
    ensure_dir(CACHE_DIR)

    # Load data
    expr_lin = load_expression(args.expr)
    meta_raw = parse_metadata(args.meta)
    meta, expr_lin = align_metadata_expression(meta_raw, expr_lin)

    # Compute entropy
    entropy = compute_entropy(expr_lin)
    anova_p, kruskal_p, variance_dict, times = test_entropy_across_time(entropy, meta)

    # Differential expression (entropy-based)
    volcano_entropy, high_cells, low_cells = differential_expression_entropy(expr_lin, entropy, top_n=10)

    # Differential expression (time-based)
    volcano_time, early_cells, late_cells = differential_expression_time(expr_lin, meta)

    # Gene name mapping
    all_genes = expr_lin.index.tolist()
    mapping_dict = map_gene_names(all_genes)
    volcano_entropy['gene_name'] = volcano_entropy['gene'].map(mapping_dict)
    volcano_time['gene_name'] = volcano_time['gene'].map(mapping_dict)

    # GO enrichment (entropy)
    logger.info("Running GO enrichment for entropy‑associated genes")
    go_entropy = run_gprofiler(
        queries=volcano_entropy.head(args.top_n)['gene'].tolist(),
        organism="scerevisiae"
    )
    if go_entropy is not None:
        go_entropy.sort_values("p_value").to_csv(out_dir / "go_entropy.csv", index=False)

    # GO enrichment (time)
    logger.info("Running GO enrichment for time‑associated genes")
    go_time = run_gprofiler(
        queries=volcano_time.head(args.top_n)['gene'].tolist(),
        organism="scerevisiae"
    )
    if go_time is not None:
        go_time.sort_values("p_value").to_csv(out_dir / "go_time.csv", index=False)

    # Cytoskeleton analysis: convert symbols to ORFs
    cyto_conv = convert_genes(CYTO_SYMBOLS)
    if cyto_conv is not None and 'converted' in cyto_conv.columns:
        cyto_orfs = cyto_conv['converted'].dropna().unique()
    else:
        logger.warning("Cytoskeleton symbol conversion failed; using empty list")
        cyto_orfs = []

    cyto_entropy, cyto_time = cytoskeleton_analysis(expr_lin, meta, cyto_orfs, mapping_dict)

    # Classification
    classification, entropy_top, time_top = classify_genes(
        volcano_entropy, volcano_time, all_genes, top_n=args.top_n
    )
    classification['gene_name'] = classification['gene'].map(mapping_dict)
    classification.to_csv(out_dir / "gene_classification.csv", index=False)

    # Cytoskeleton classification
    cyto_class = classification[classification['gene'].isin(cyto_orfs)]
    cyto_class.to_csv(out_dir / "cyto_classification.csv", index=False)

    # Save CSVs
    volcano_entropy.to_csv(out_dir / "entropy_volcano.csv", index=False)
    volcano_time.to_csv(out_dir / "time_volcano.csv", index=False)
    cyto_time.to_csv(out_dir / "cyto_time_expression.csv", index=False)
    pd.DataFrame(list(variance_dict.items()), columns=['time','variance']).to_csv(out_dir / "entropy_variance.csv", index=False)

    # Generate plots
    entropy_box_path = out_dir / "entropy_box.png"
    plot_entropy_box(entropy, meta, entropy_box_path)

    cyto_entropy_plot_path = out_dir / "cyto_entropy.png"
    plot_cyto_entropy(cyto_entropy, meta, cyto_entropy_plot_path)

    variance_plot_path = out_dir / "variance.png"
    plot_variance_trend(variance_dict, times, variance_plot_path)

    volcano_entropy_path = out_dir / "entropy_volcano.png"
    plot_volcano(volcano_entropy, highlight_genes=cyto_orfs,
                 title="High‑ vs Low‑Entropy Cells", out_path=volcano_entropy_path)

    volcano_time_path = out_dir / "time_volcano.png"
    plot_volcano(volcano_time, highlight_genes=cyto_orfs,
                 title="Late vs Early Time Points", out_path=volcano_time_path)

    cyto_traj_path = out_dir / "cyto_time_trends.png"
    plot_cyto_trajectories(cyto_time, cyto_traj_path)

    # Prepare results dictionary for PDF
    results_for_pdf = {
        'anova_p': anova_p,
        'kruskal_p': kruskal_p,
        'entropy_box': entropy_box_path,
        'cyto_entropy_plot': cyto_entropy_plot_path,
        'variance_plot': variance_plot_path,
        'volcano_entropy': volcano_entropy_path,
        'volcano_time': volcano_time_path,
        'volcano_entropy_df': volcano_entropy,
        'volcano_time_df': volcano_time,
        'go_entropy': go_entropy,
        'go_time': go_time,
        'cyto_trajectories': cyto_traj_path,
        'classification': classification,
        'entropy_top': entropy_top,
        'time_top': time_top
    }

    # Generate PDF report
    pdf_path = out_dir / "entropy_report.pdf"
    generate_pdf_report(pdf_path, results_for_pdf)

    # Terminal summary
    logger.info("="*50)
    logger.info("ANALYSIS COMPLETE")
    logger.info(f"Outputs saved in {out_dir.absolute()}")
    logger.info(f"Global entropy ANOVA p = {anova_p:.3e}")
    logger.info(f"Number of significant genes (entropy DE, p<{PVAL_THRESH}): {(volcano_entropy['pval'] < PVAL_THRESH).sum()}")
    logger.info(f"Number of significant genes (time DE, p<{PVAL_THRESH}): {(volcano_time['pval'] < PVAL_THRESH).sum()}")
    logger.info(f"Overlap of top {args.top_n} genes: {len(entropy_top.intersection(time_top))}")
    logger.info("Cytoskeleton gene classification:")
    logger.info(cyto_class['class'].value_counts().to_string())
    logger.info("="*50)

if __name__ == "__main__":
    main()