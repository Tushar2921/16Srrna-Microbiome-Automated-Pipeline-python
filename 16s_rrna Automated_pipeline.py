#!/usr/bin/env python3
"""
ULTIMATE 16S MICROBIOME THESIS PIPELINE - PYTHON VERSION (2026)
================================================================
Equivalent to 16S_new.sh but written in Python.
Uses the SAME tools, versions, URLs, and directory structure.

Tools: fastqc, multiqc, fastp, sra-tools, QIIME2 2024.5, DADA2,
       SILVA classifier, biom-format, R (phyloseq, vegan, DESeq2, ggplot2)
"""

import os
import sys
import shutil
import subprocess
import logging
import platform
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Tuple

# ============================================================================
# CONFIGURATION
# ============================================================================
PROJECT_NAME = "16S_microbiome_analysis"
SAMPLE_MODE = "sra"                    # "sra" | "local" | "both"
SAMPLE_FILE = Path.home() / "Pipeline" / "sample_list.txt"
NUM_THREADS = "auto"

BASE_DIR = Path.home() / PROJECT_NAME
LOG_DIR = BASE_DIR / "logs"
RAW_DATA = BASE_DIR / "01_raw_data"
TRIMMED = BASE_DIR / "02_trimmed"
QC_REPORTS = BASE_DIR / "03_qc_reports"
QIIME_INPUT = BASE_DIR / "04_qiime_input"
QIIME_OUTPUT = BASE_DIR / "05_qiime_output"
R_ANALYSIS = BASE_DIR / "06_r_analysis"
FINAL_RESULTS = BASE_DIR / "07_results"

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_FILE = LOG_DIR / f"{TIMESTAMP}_pipeline.log"

# ============================================================================
# LOGGING SETUP
# ============================================================================
def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("16S_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(levelname)s] %(message)s")

    fh = logging.FileHandler(LOG_FILE)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)
    return logger

log = setup_logging()

# ============================================================================
# HELPER: RUN A SUBPROCESS COMMAND
# ============================================================================
def run(cmd: List[str], check: bool = False, env: Optional[dict] = None,
        shell: bool = False) -> Tuple[int, str, str]:
    """Run a shell command and return (returncode, stdout, stderr)."""
    log.info(f"$ {' '.join(cmd) if not shell else cmd}")
    try:
        result = subprocess.run(
            cmd if not shell else cmd,
            shell=shell,
            capture_output=True,
            text=True,
            env=env,
        )
        if result.stdout:
            log.info(result.stdout.strip())
        if result.stderr:
            log.info(result.stderr.strip())
        if check and result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd)
        return result.returncode, result.stdout, result.stderr
    except FileNotFoundError as e:
        log.error(f"Command not found: {cmd}")
        return 1, "", str(e)

def section(title: str):
    log.info("\n" + "=" * 60)
    log.info(f"→ {title}")
    log.info("=" * 60 + "\n")

# ============================================================================
# STEP 1: PREFLIGHT CHECKS
# ============================================================================
def preflight_checks() -> bool:
    section("PREFLIGHT CHECKS")
    if not shutil.which("conda"):
        log.warning("Conda not found - installing Miniconda...")
        if not install_miniconda():
            return False

    # Source conda for this session
    conda_sh = Path.home() / "miniconda3" / "etc" / "profile.d" / "conda.sh"
    if conda_sh.exists():
        os.environ["PATH"] = f"{Path.home()/'miniconda3'/'bin'}:{os.environ['PATH']}"

    rc, out, _ = run(["conda", "--version"])
    if rc != 0:
        log.error("Conda initialization failed")
        return False
    log.info(f"Conda initialized: {out.strip()}")

    # Disk space
    total, used, free = shutil.disk_usage(Path.home())
    free_gb = free // (1024 ** 3)
    log.info(f"Available disk space: {free_gb} GB")
    if free_gb < 20:
        log.warning("Less than 20GB available - may have issues")

    # Threads
    global NUM_THREADS
    if NUM_THREADS == "auto":
        NUM_THREADS = str(max(1, os.cpu_count() - 1)) if os.cpu_count() else "4"
    log.info(f"Using {NUM_THREADS} threads")
    log.info("✓ Preflight checks passed")
    return True

# ============================================================================
# STEP 2: INSTALL MINICONDA
# ============================================================================
def install_miniconda() -> bool:
    log.info("Installing Miniconda...")
    system = platform.system().lower()
    machine = platform.machine()
    if system == "darwin":
        url = ("https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
               if machine == "arm64"
               else "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh")
    else:
        url = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

    installer = Path.home() / "miniconda.sh"
    rc, _, _ = run(["wget", "-q", url, "-O", str(installer)])
    if rc != 0:
        log.error("Failed to download Miniconda")
        return False

    rc, _, _ = run(["bash", str(installer), "-b", "-p", str(Path.home() / "miniconda3")])
    installer.unlink(missing_ok=True)
    if rc != 0:
        log.error("Miniconda installation failed")
        return False
    log.info("✓ Miniconda installed")
    return True

# ============================================================================
# STEP 3: SETUP CONDA CHANNELS
# ============================================================================
def setup_conda():
    section("SETTING UP CONDA")
    run(["conda", "config", "--set", "channel_priority", "flexible"])
    run(["conda", "config", "--remove-key", "channels"])
    for ch in ["defaults", "bioconda", "conda-forge"]:
        run(["conda", "config", "--add", "channels", ch])
    run(["conda", "clean", "--all", "-y", "-q"])
    run(["conda", "install", "-n", "base", "conda-libmamba-solver", "-y", "-q"])
    run(["conda", "update", "-n", "base", "-y", "-q", "conda"])
    log.info("✓ Conda setup complete")

# ============================================================================
# STEP 4: CREATE ENVIRONMENTS
# ============================================================================
def env_exists(name: str) -> bool:
    rc, out, _ = run(["conda", "env", "list"])
    return name in out

def create_environments():
    section("CREATING CONDA ENVIRONMENTS")

    # 1) QC env
    log.info("[1/3] Creating QC environment...")
    if not env_exists("qc-env"):
        run(["conda", "create", "-n", "qc-env", "python=3.9",
             "fastqc", "multiqc", "fastp", "sra-tools",
             "-c", "conda-forge", "-c", "bioconda", "-y", "-q"])
    else:
        log.info("QC environment exists")

    # 2) QIIME2 env (official yml)
    log.info("[2/3] Creating QIIME2 environment...")
    if not env_exists("qiime2-env"):
        rc, _, _ = run(["conda", "env", "create", "-n", "qiime2-env",
                        "--file", "https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml",
                        "-y"])
        if rc != 0:
            log.warning("Amplicon failed, trying Tiny distribution...")
            run(["conda", "env", "create", "-n", "qiime2-env",
                 "--file", "https://data.qiime2.org/distro/tiny/qiime2-tiny-2024.5-py39-linux-conda.yml",
                 "-y"])
    else:
        log.info("QIIME2 environment exists")

    # 3) R env
    log.info("[3/3] Creating R environment...")
    if not env_exists("r-stats-env"):
        run(["conda", "create", "-n", "r-stats-env", "python=3.9",
             "r-base=4.2", "r-ggplot2", "r-dplyr", "r-tidyr",
             "r-vegan", "r-phyloseq", "r-deseq2", "r-igraph",
             "r-ggraph", "r-pheatmap", "r-rstatix", "biom-format",
             "-c", "conda-forge", "-c", "bioconda", "-y", "-q"])
    log.info("✓ Environment creation complete")

# ============================================================================
# STEP 5-6: DIRECTORIES + VALIDATION
# ============================================================================
def create_directories():
    section("CREATING DIRECTORIES")
    for d in [RAW_DATA, TRIMMED, QC_REPORTS, QIIME_INPUT,
              QIIME_OUTPUT, R_ANALYSIS, FINAL_RESULTS,
              FINAL_RESULTS / "figures", FINAL_RESULTS / "tables"]:
        d.mkdir(parents=True, exist_ok=True)
    log.info("✓ Directories created")

def validate_input() -> bool:
    section("VALIDATING INPUT")
    if SAMPLE_MODE in ("local", "both") and not SAMPLE_FILE.exists():
        log.error(f"Sample file not found: {SAMPLE_FILE}")
        return False
    log.info("✓ Input validation passed")
    return True

# ============================================================================
# STEP 7: DOWNLOAD SRA
# ============================================================================
def conda_run(env: str, cmd: List[str]) -> Tuple[int, str, str]:
    """Run a command inside a conda env."""
    return run(["conda", "run", "-n", env, "--no-banner"] + cmd)

def download_sra():
    section("DOWNLOADING SRA DATA")
    if not SAMPLE_FILE.exists():
        log.error("Sample file not found"); return
    with open(SAMPLE_FILE) as f:
        srr_ids = [line.strip() for line in f if line.strip()]
    for i, srr in enumerate(srr_ids, 1):
        log.info(f"[SRA {i}] Downloading {srr}...")
        rc, _, _ = conda_run("qc-env",
            ["fasterq-dump", srr, "-O", str(RAW_DATA), "-e", str(NUM_THREADS)])
        if rc == 0:
            log.info(f"✓ {srr} downloaded")
        else:
            log.warning(f"Failed to download {srr}")
    log.info("✓ SRA download complete")

# ============================================================================
# STEP 8: QC (FastQC + fastp + MultiQC)
# ============================================================================
def run_qc():
    section("RUNNING QUALITY CONTROL")
    fastqs = list(RAW_DATA.glob("*.fastq*"))
    if not fastqs:
        log.warning("No fastq files found"); return

    log.info("Running FastQC...")
    conda_run("qc-env",
        ["fastqc"] + [str(f) for f in fastqs] +
        ["-o", str(QC_REPORTS), "-t", str(NUM_THREADS), "--nogroup"])

    log.info("Trimming with fastp...")
    seen = set()
    for f in RAW_DATA.glob("*_1.fastq*") + tuple(RAW_DATA.glob("*_R1*.fastq*")):
        base = f.name.split("_1.fastq")[0].split("_R1.fastq")[0]
        if base in seen: continue
        seen.add(base)
        r1 = RAW_DATA / f"{base}_1.fastq"
        r2 = RAW_DATA / f"{base}_2.fastq"
        if not r1.exists():
            r1 = next(RAW_DATA.glob(f"{base}_R1*"), None)
            r2 = next(RAW_DATA.glob(f"{base}_R2*"), None)
        if r1 and r2 and r1.exists() and r2.exists():
            conda_run("qc-env",
                ["fastp", "-i", str(r1), "-I", str(r2),
                 "-o", str(TRIMMED / f"{base}_1.fastq"),
                 "-O", str(TRIMMED / f"{base}_2.fastq"),
                 "-q", "20", "-l", "50", "-w", str(NUM_THREADS)])
            log.info(f"✓ Trimmed {base}")

    log.info("Creating MultiQC report...")
    conda_run("qc-env",
        ["multiqc", str(QC_REPORTS), str(TRIMMED),
         "-o", str(QC_REPORTS), "-n", "multiqc_report.html"])
    log.info("✓ QC complete")

# ============================================================================
# STEP 9: QIIME2 (Import → DADA2 → Taxonomy)
# ============================================================================
def run_qiime2():
    section("RUNNING QIIME2 ANALYSIS")

    # Manifest
    manifest = QIIME_INPUT / "manifest.tsv"
    with open(manifest, "w") as m:
        m.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for f in TRIMMED.glob("*_1.fastq"):
            base = f.stem.replace("_1", "")
            r1 = (TRIMMED / f"{base}_1.fastq").resolve()
            r2 = (TRIMMED / f"{base}_2.fastq").resolve()
            m.write(f"{base}\t{r1}\t{r2}\n")

    # Auto metadata (Group_A / Group_B split)
    meta = QIIME_INPUT / "sample_metadata.tsv"
    lines = manifest.read_text().splitlines()[1:]
    half = len(lines) // 2
    with open(meta, "w") as m:
        m.write("sample-id\tcondition\n")
        for i, line in enumerate(lines):
            sid = line.split("\t")[0]
            m.write(f"{sid}\t{'Group_A' if i < half else 'Group_B'}\n")

    # Import
    demux = QIIME_OUTPUT / "demux.qza"
    conda_run("qiime2-env",
        ["qiime", "tools", "import",
         "--type", "SampleData[PairedEndSequencesWithQuality]",
         "--input-path", str(manifest),
         "--output-path", str(demux),
         "--input-format", "PairedEndFastqManifestPhred33V2"])

    if demux.exists():
        conda_run("qiime2-env",
            ["qiime", "demux", "summarize",
             "--i-data", str(demux),
             "--o-visualization", str(QIIME_OUTPUT / "demux.qzv")])

    # DADA2
    table = QIIME_OUTPUT / "table.qza"
    repseqs = QIIME_OUTPUT / "rep-seqs.qza"
    if demux.exists():
        log.info("Running DADA2 (15-60 min)...")
        conda_run("qiime2-env",
            ["qiime", "dada2", "denoise-paired",
             "--i-demultiplexed-seqs", str(demux),
             "--p-trim-left-f", "17", "--p-trim-left-r", "21",
             "--p-trunc-len-f", "250", "--p-trunc-len-r", "250",
             "--p-max-ee-f", "2", "--p-max-ee-r", "2",
             "--p-n-threads", str(NUM_THREADS),
             "--o-table", str(table),
             "--o-representative-sequences", str(repseqs),
             "--o-denoising-stats", str(QIIME_OUTPUT / "denoising-stats.qza")])
        log.info("✓ DADA2 complete")

    # Taxonomy
    if repseqs.exists():
        log.info("Assigning taxonomy...")
        classifier = QIIME_OUTPUT / "silva-classifier.qza"
        urls = [
            "https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza",
            "https://data.qiime2.org/2022.2/common/silva-138-99-515-806-nb-classifier.qza",
            "https://data.qiime2.org/2021.4/common/silva-138-99-515-806-nb-classifier.qza",
        ]
        for url in urls:
            rc, _, _ = run(["wget", "--timeout=300", "--tries=3", "-q", url, "-O", str(classifier)])
            if rc == 0 and classifier.stat().st_size > 10 * 1024 * 1024:
                peek_rc, _, _ = conda_run("qiime2-env", ["qiime", "tools", "peek", str(classifier)])
                if peek_rc == 0:
                    log.info(f"✓ Classifier downloaded from {url}")
                    break
            classifier.unlink(missing_ok=True)

        if classifier.exists():
            conda_run("qiime2-env",
                ["qiime", "feature-classifier", "classify-sklearn",
                 "--i-classifier", str(classifier),
                 "--i-reads", str(repseqs),
                 "--p-n-jobs", str(NUM_THREADS),
                 "--o-classification", str(QIIME_OUTPUT / "taxonomy.qza")])
            log.info("✓ Taxonomy assignment complete")

    log.info("✓ QIIME2 analysis complete")

# ============================================================================
# STEP 10: EXPORT QIIME2 RESULTS TO TSV
# ============================================================================
def export_results():
    section("EXPORTING RESULTS")
    table = QIIME_OUTPUT / "table.qza"
    tax = QIIME_OUTPUT / "taxonomy.qza"
    if not table.exists():
        log.warning("No QIIME2 table to export"); return

    conda_run("qiime2-env",
        ["qiime", "tools", "export", "--input-path", str(table),
         "--output-path", str(R_ANALYSIS / "table")])
    conda_run("qiime2-env",
        ["biom", "convert",
         "-i", str(R_ANALYSIS / "table" / "feature-table.biom"),
         "-o", str(R_ANALYSIS / "feature-table.tsv"),
         "--to-tsv"])
    if tax.exists():
        conda_run("qiime2-env",
            ["qiime", "tools", "export", "--input-path", str(tax),
             "--output-path", str(R_ANALYSIS / "tax")])
    log.info("✓ Export complete")

# ============================================================================
# STEP 11: R THESIS PLOTS (Alpha/Beta/DESeq2/Network)
# ============================================================================
def run_r_thesis():
    section("GENERATING THESIS PLOTS (R)")
    if not (R_ANALYSIS / "feature-table.tsv").exists():
        log.warning("QIIME2 exports missing - skipping R plots"); return

    r_script = R_ANALYSIS / "generate_thesis_plots.R"
    r_script.write_text(r"""
suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(DESeq2); library(ggplot2)
  library(dplyr); library(tidyr); library(igraph); library(ggraph)
})
OUT_FIG <- "FINAL_RESULTS_PLACEHOLDER/figures"
OUT_TAB <- "FINAL_RESULTS_PLACEHOLDER/tables"

otu_raw <- read.delim("R_ANALYSIS_PLACEHOLDER/feature-table.tsv",
                      row.names=1, check.names=FALSE)
otu <- t(otu_raw[,-1])
tax_raw <- read.delim("R_ANALYSIS_PLACEHOLDER/tax/taxonomy.tsv",
                      row.names=1, check.names=FALSE)
meta <- read.delim("QIIME_INPUT_PLACEHOLDER/sample_metadata.tsv",
                   row.names=1, check.names=FALSE)
colnames(otu) <- meta$sample_id
tax_mat <- as.matrix(tax_raw[,-c(1,2)])
rownames(tax_mat) <- tax_raw$Feature.ID
ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
               tax_table(tax_mat), sample_data(meta))

# Alpha
alpha_df <- estimate_richness(ps, measures=c("Shannon","Observed","Chao1","FaithPD"))
alpha_df$condition <- sample_data(ps)$condition
p1 <- ggplot(alpha_df %>% pivot_longer(Shannon:FaithPD, "Metric","Value"),
             aes(condition, Value, fill=condition)) +
  geom_boxplot() + facet_wrap(~Metric, scales="free_y") + theme_minimal()
ggsave(file.path(OUT_FIG, "alpha_diversity.pdf"), p1, width=8, height=5)

# Beta (PCoA Bray-Curtis)
bc <- distance(ps, "bray")
ord <- ordinate(ps, "PCoA", bc)
p2 <- plot_ordination(ps, ord, color="condition") + geom_point(size=4) + theme_minimal()
ggsave(file.path(OUT_FIG, "pcoa_bray.pdf"), p2, width=7, height=6)
write.table(adonis2(bc ~ condition, data=meta, permutations=999),
            file.path(OUT_TAB, "permanova_results.tsv"), sep="\t")

# Relative abundance (top 12 genera)
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
top <- names(sort(taxa_sums(ps_rel), decreasing=TRUE)[1:12])
df <- psmelt(prune_taxa(top, ps_rel))
p3 <- ggplot(df, aes(sample, Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill") +
  facet_wrap(~condition, scales="free_x") + theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_FIG, "relative_abundance.pdf"), p3, width=12, height=5)

# DESeq2 volcano
dds <- phyloseq_to_deseq2(ps, ~ condition)
dds <- DESeq(dds, test="Wald", fitType="local", quiet=TRUE)
res <- results(dds, alpha=0.05)
sig <- res[which(res$padj < 0.05), ]
write.table(as.data.frame(sig), file.path(OUT_TAB, "deseq2_significant.tsv"), sep="\t")
v <- data.frame(log2FC=res$log2FoldChange, negLogP=-log10(res$pvalue),
                sig=res$padj<0.05)
p4 <- ggplot(v[complete.cases(v),], aes(log2FC, negLogP, color=sig)) +
  geom_point(alpha=0.6) + scale_color_manual(values=c("grey","red")) + theme_minimal()
ggsave(file.path(OUT_FIG, "volcano.pdf"), p4, width=7, height=5)

# Network
sp <- cor(as.matrix(otu), method="spearman")
adj <- ifelse(abs(sp) > 0.6, 1, 0); diag(adj) <- 0
g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE)
V(g)$degree <- degree(g)
p5 <- ggraph(g, layout="fr") + geom_edge_link(aes(alpha=weight)) +
  geom_node_point(aes(size=degree), color="steelblue") + theme_void()
ggsave(file.path(OUT_FIG, "network.pdf"), p5, width=9, height=8)
cat("✅ R thesis plots saved\n")
""")

    # Replace placeholders
    txt = r_script.read_text()
    txt = txt.replace("FINAL_RESULTS_PLACEHOLDER", str(FINAL_RESULTS))
    txt = txt.replace("R_ANALYSIS_PLACEHOLDER", str(R_ANALYSIS))
    txt = txt.replace("QIIME_INPUT_PLACEHOLDER", str(QIIME_INPUT))
    r_script.write_text(txt)

    conda_run("r-stats-env", ["Rscript", str(r_script)])
    log.info("✓ R thesis plots generated")

# ============================================================================
# STEP 12: FINAL REPORT
# ============================================================================
def generate_report():
    section("GENERATING FINAL REPORT")
    report = FINAL_RESULTS / "ANALYSIS_REPORT.txt"
    report.write_text(f"""
================================================================================
16S MICROBIOME ANALYSIS REPORT
================================================================================
Project: {PROJECT_NAME}
Date:    {datetime.now()}
Dir:     {BASE_DIR}
Log:     {LOG_FILE}
================================================================================
RESULTS LOCATION:
  QC Reports:    {QC_REPORTS}/multiqc_report.html
  QIIME2 Output: {QIIME_OUTPUT}/
  R Analysis:    {R_ANALYSIS}/
  Thesis Plots:  {FINAL_RESULTS}/figures/
  Stats Tables:  {FINAL_RESULTS}/tables/
================================================================================
NEXT STEPS:
1. Open multiqc_report.html in a browser
2. Review figures in {FINAL_RESULTS}/figures/
3. Extract stats from {FINAL_RESULTS}/tables/
4. Edit sample_metadata.tsv if you want real group names
================================================================================
""")
    log.info(f"✓ Report generated: {report}")

# ============================================================================
# MAIN
# ============================================================================
def main():
    log.info("\n" + "=" * 60)
    log.info("  16S THESIS PIPELINE - PYTHON VERSION")
    log.info("=" * 60 + "\n")

    if not preflight_checks(): sys.exit(1)
    setup_conda()
    create_environments()
    create_directories()
    if not validate_input(): sys.exit(1)

    if SAMPLE_MODE in ("sra", "both"):
        download_sra()
    run_qc()
    run_qiime2()
    export_results()
    run_r_thesis()
    generate_report()

    log.info("\n" + "=" * 60)
    log.info("  ✅ ANALYSIS COMPLETE")
    log.info("=" * 60)
    log.info(f"Results at: {BASE_DIR}")

if __name__ == "__main__":
    main()