# Yeast Aging Transcriptomics Analysis

## Overview
This project explores transcriptional changes associated with yeast aging using bulk RNA-seq data. The analysis focuses on identifying patterns in gene expression dynamics, entropy-based variability, and functional enrichment across time.

The pipeline integrates data preprocessing, statistical analysis, entropy calculations, and visualization to derive biological insights into aging-related processes.

---

## Dataset
- Source: GEO (Gene Expression Omnibus)
- Example dataset: GSE210032
- Data includes:
  - Normalized gene expression matrices
  - Metadata across time points

---

## Objectives
- Analyze temporal gene expression changes during aging
- Compute entropy as a measure of transcriptional variability
- Identify genes with significant dynamic behavior
- Perform GO enrichment analysis
- Visualize trends and distributions

---

## Pipeline Overview

1. **Data Acquisition**
   - Download GEO datasets
   - Extract and organize expression matrices

2. **Preprocessing**
   - Normalization (if required)
   - Gene filtering

3. **Entropy Analysis**
   - Compute entropy for each gene across time
   - Identify high-variability genes

4. **Statistical Analysis**
   - Variance and differential trends
   - Volcano plots

5. **Functional Enrichment**
   - GO enrichment using g:Profiler

6. **Visualization**
   - Time-series plots
   - Entropy distributions
   - Volcano plots

---

## Project Structure

```
yeast_aging/
│── data/                  # Raw and processed datasets
│── outputs/               # Analysis results (plots, tables)
│── outputsGSE210032/      # Dataset-specific outputs
│── yeast_aging.py         # Main analysis script
│── analysis.log           # Logging output
```

---

## Key Outputs

- `entropy_volcano.png` – Genes with significant entropy changes
- `cyto_time_trends.png` – Temporal expression trends
- `go_entropy.csv` – GO enrichment results
- `entropy_report.pdf` – Summary report

---

## Requirements

- Python 3.x
- Libraries:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - gprofiler-official

Install dependencies:

```bash
pip install pandas numpy matplotlib seaborn gprofiler-official
```

---

## Usage

Run the main script:

```bash
python yeast_aging.py
```

Outputs will be generated in the `outputs/` directory.

---

## Notes

- Large files (logs, cache, intermediate outputs) should be excluded using `.gitignore`
- Ensure consistent gene identifiers for enrichment analysis

---

## Future Directions

- Integrate single-cell RNA-seq analysis
- Apply trajectory inference methods
- Compare across different yeast strains
- Incorporate evolutionary conservation analysis

---

## Author
Anirudh Baliga

---

## License
This project is for academic and research purposes.
