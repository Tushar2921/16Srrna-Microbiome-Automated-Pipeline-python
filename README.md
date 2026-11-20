# 🧬 16S Microbiome Analysis Pipeline — Python Version

> **A complete, automated, thesis-ready 16S rRNA amplicon sequencing pipeline** written in Python. Orchestrates bioinformatics tools through `subprocess`, providing clean logging, modular design, and cross-platform compatibility.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Status](https://img.shields.io/badge/Status-Production%20Ready-green.svg)]()

---

## 📋 Table of Contents
1. [Overview](#-overview)
2. [Why Python?](#-why-python)
3. [Features](#-features)
4. [Prerequisites](#-prerequisites)
5. [Installation & Setup](#-installation--setup)
6. [How to Run](#-how-to-run)
7. [Configuration](#-configuration)
8. [Output Structure](#-output-structure)
9. [Understanding Each Output](#-understanding-each-output)
10. [Further Analysis](#-further-analysis)
11. [Troubleshooting](#-troubleshooting)
12. [Citation](#-citation)

---

## 🎯 Overview

This Python pipeline automates the **entire 16S microbiome workflow**:

```
Raw SRA Data → QC → Trimming → DADA2 → Taxonomy → Diversity → 
Differential Abundance → Networks → Publication Figures → Thesis Report
```

It's the **Python equivalent** of the Bash version, with:
- Better structured logging
- Cleaner error handling
- Easier to extend and customize
- Cross-platform support (Linux, macOS, Windows WSL)

---

## 🐍 Why Python?

| Feature | Bash Version | Python Version |
|---------|--------------|----------------|
| **Logging** | Basic `tee` | Structured `logging` module |
| **Error Handling** | Exit codes | Exceptions + try/except |
| **Extensibility** | Harder | Easy (import & extend) |
| **Cross-Platform** | Linux/macOS | Linux/macOS/Windows (WSL) |
| **Readability** | Shell syntax | Clean Python code |

---

## ✨ Features

- ✅ **Python-native** — clean, maintainable, well-documented
- ✅ **Structured logging** — timestamped logs with severity levels
- ✅ **Zero manual steps** — fully automated end-to-end
- ✅ **Modular design** — each step is a separate function
- ✅ **Publication-ready outputs** — PDFs, TSVs, statistical tables
- ✅ **Type hints** — better IDE support and code clarity

---

## 🛠️ Prerequisites

| Requirement | Details |
|-------------|---------|
| **Python** | 3.9 or higher |
| **OS** | Linux, macOS, or Windows (WSL) |
| **Disk Space** | Minimum 20 GB free |
| **RAM** | Minimum 8 GB (16 GB recommended) |
| **Internet** | Required for downloading tools & classifiers |

---

## 🚀 Installation & Setup

### Step 1: Download the Pipeline

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/16S-Microbiome-Pipeline-Python.git
cd 16S-Microbiome-Pipeline-Python

# OR: Download manually from GitHub
# Click "Code" → "Download ZIP" → Extract
```

### Step 2: Install Python Dependencies

```bash
pip install -r requirements.txt
```

> 💡 The `requirements.txt` is minimal — most tools are installed via conda automatically.

### Step 3: Prepare Your Sample List

Create `~/Pipeline/sample_list.txt` with SRA accessions (one per line):

```bash
mkdir -p ~/Pipeline
nano ~/Pipeline/sample_list.txt
```

**Example content:**
```
SRR12345678
SRR12345679
SRR12345680
```

### Step 4: Make Executable (Optional)

```bash
chmod +x 16s_microbiome_pipeline.py
```

---

## ▶️ How to Run

```bash
python3 16s_microbiome_pipeline.py
```

**That's it!** The pipeline will:
1. Auto-install Miniconda (if missing)
2. Create 3 conda environments (`qc-env`, `qiime2-env`, `r-stats-env`)
3. Download your SRA data
4. Run QC, trimming, DADA2, taxonomy
5. Generate all figures and statistics
6. Output a thesis interpretation guide

**Total runtime:** 2–6 hours (depends on sample count & hardware)

---

## ⚙️ Configuration

Edit these variables at the **top of the script**:

```python
PROJECT_NAME = "16S_microbiome_analysis"   # Output folder name
SAMPLE_MODE = "sra"                         # "sra" | "local" | "both"
SAMPLE_FILE = Path.home() / "Pipeline" / "sample_list.txt"
NUM_THREADS = "auto"                        # Auto-detects CPU cores
```

### For Local FASTQ Files (instead of SRA)

```python
SAMPLE_MODE = "local"
```
Then place your FASTQ files in `~/16S_microbiome_analysis/01_raw_data/`

---

## 📂 Output Structure

After running, you'll find:

```
~/16S_microbiome_analysis/
│
├── 01_raw_data/              # Downloaded SRA files
├── 02_trimmed/               # Quality-filtered reads
├── 03_qc_reports/            # FastQC + MultiQC reports
│   └── multiqc_report.html   ← Open in browser
├── 04_qiime_input/           # Manifest & metadata
├── 05_qiime_output/          # QIIME2 artifacts (.qza/.qzv)
├── 06_r_analysis/            # Exported tables for R
├── 07_results/               # 🎯 FINAL THESIS-READY OUTPUTS
│   ├── figures/              # Publication PDFs
│   │   ├── alpha_diversity.pdf
│   │   ├── pcoa_bray.pdf
│   │   ├── relative_abundance.pdf
│   │   ├── volcano.pdf
│   │   └── network.pdf
│   ├── tables/               # Statistical tables (TSV)
│   │   ├── alpha_kruskal_pval.txt
│   │   ├── permanova_results.tsv
│   │   └── deseq2_significant.tsv
│   └── THESIS_READY_GUIDE.txt
└── logs/                     # Pipeline execution logs
```

---

## 🔍 Understanding Each Output

### 📊 1. Quality Control (`03_qc_reports/multiqc_report.html`)

**What it is:** Aggregated QC report showing read quality, GC content, adapter contamination.

**How to use:**
```bash
firefox 03_qc_reports/multiqc_report.html
```

**Thesis phrasing:**
> "Raw reads were quality-filtered using fastp (Q≥20, min length 50 bp). MultiQC confirmed >85% read retention with Q30 scores met across all samples."

---

### 📈 2. Alpha Diversity (`07_results/figures/alpha_diversity.pdf`)

**What it is:** Boxplots of Shannon, Chao1, Observed ASVs, FaithPD by group.

**How to interpret:**
- Higher values = more diverse community
- Compare Group_A vs Group_B
- Check `tables/alpha_kruskal_pval.txt` for significance

**Thesis phrasing:**
> "Alpha diversity (Shannon index) did not differ significantly between groups (Kruskal-Wallis, p=0.XX), indicating comparable within-sample complexity."

---

### 🌍 3. Beta Diversity / PCoA (`07_results/figures/pcoa_bray.pdf`)

**What it is:** Principal Coordinates Analysis showing community clustering.

**How to interpret:**
- Points close together = similar communities
- Check `tables/permanova_results.tsv` for R² and p-value

**Thesis phrasing:**
> "PCoA ordination revealed distinct clustering by experimental condition (PERMANOVA, R²=0.XX, p<0.01)."

---

### 🧬 4. Relative Abundance (`07_results/figures/relative_abundance.pdf`)

**What it is:** Stacked barplot of top 12 genera per sample.

**How to interpret:**
- Identify dominant taxa
- Spot group-specific shifts

**Thesis phrasing:**
> "The microbiome was dominated by [Genus1] (X%), [Genus2] (Y%). Group-specific shifts were observed with [Taxon] enriched in Group_A."

---

### 🔥 5. Volcano Plot (`07_results/figures/volcano.pdf`)

**What it is:** DESeq2 differential abundance results.

**How to interpret:**
- X-axis: log2 fold change (enrichment/depletion)
- Y-axis: -log10(p-value)
- Red points: significant taxa (padj < 0.05)

**Thesis phrasing:**
> "DESeq2 identified X taxa significantly differentially abundant (padj < 0.05), including [Taxon1] (log2FC=2.1)."

---

### 🕸️ 6. Network (`07_results/figures/network.pdf`)

**What it is:** Spearman co-occurrence network (|r| > 0.6).

**How to interpret:**
- Large nodes = keystone taxa (high degree)
- Edges = correlations (positive/negative)

**Thesis phrasing:**
> "Co-occurrence network revealed modular structure. [Taxon] exhibited highest degree centrality, suggesting a keystone ecological role."

---

## 🔬 Further Analysis

### Functional Prediction (PICRUSt2)

```bash
conda activate picrust-env
picrust2_pipeline.py \
  -s 05_qiime_output/rep-seqs.fasta \
  -i 06_r_analysis/feature-table.biom \
  -o picrust2_output \
  -p 8
```

### Phylogenetic Tree (for UniFrac)

```bash
conda activate qiime2-env
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 05_qiime_output/rep-seqs.qza \
  --o-alignment aligned.qza \
  --o-masked-alignment masked.qza \
  --o-tree unrooted.qza \
  --o-rooted-tree rooted.qza
```

### ANCOM (Alternative Differential Abundance)

```bash
qiime composition ancom \
  --i-table 05_qiime_output/table.qza \
  --m-metadata-file 04_qiime_input/sample_metadata.tsv \
  --m-metadata-column condition \
  --output-dir ancom_results
```

---

## 🐛 Troubleshooting

### ❌ Classifier Download Fails

**Symptom:** `Could not download classifier automatically`

**Fix:**
```bash
wget -O ~/16S_microbiome_analysis/05_qiime_output/silva-classifier.qza \
  https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza
```
Then re-run the pipeline.

---

### ❌ Conda Environment Issues

**Fix:**
```bash
conda env remove -n qiime2-env
conda env remove -n qc-env
conda env remove -n r-stats-env
python3 16s_microbiome_pipeline.py
```

---

### ❌ DADA2 Fails with "No reads retained"

**Cause:** Truncation parameters too aggressive.

**Fix:** Edit the script — reduce `--p-trunc-len-f` and `--p-trunc-len-r` to 200.

---

### ❌ R Plots Not Generated

**Check:** `logs/*_pipeline.log` for R errors.

**Fix:** Ensure `r-stats-env` was created successfully:
```bash
conda activate r-stats-env
Rscript --version
```

---

### ❌ "No fastq files found"

**Cause:** SRA download failed or sample list is empty.

**Fix:**
1. Check `~/Pipeline/sample_list.txt` has valid SRR IDs
2. Verify internet connection
3. Try manual download: `fasterq-dump SRR12345678 -O ~/16S_microbiome_analysis/01_raw_data/`

---

### ❌ Python Version Errors

**Symptom:** `SyntaxError` or `ImportError`

**Fix:** Ensure Python 3.9+:
```bash
python3 --version
# Should show 3.9.x or higher
```

---

## 📚 Citation

If you use this pipeline, please cite the underlying tools:

```bibtex
@article{bolyen2019qiime2,
  title={Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2},
  author={Bolyen, Evan and others},
  journal={Nature biotechnology},
  year={2019}
}

@article{callahan2016dada2,
  title={DADA2: High-resolution sample resolution from Illumina amplicon data},
  author={Callahan, Benjamin J and others},
  journal={Nature methods},
  year={2016}
}
```

---

## 🤝 Contributing

Contributions welcome! Open an issue or submit a pull request.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) file.

---

## 👤 Author

**Tushar**  
Bioinformatics | Clinical Genomics | MEDBIORA  
🔗 [LinkedIn](YOUR_LINKEDIN_URL) | ✉️ [Email](mailto:YOUR_EMAIL)

---

⭐ **Star this repo if it helped your research!**
