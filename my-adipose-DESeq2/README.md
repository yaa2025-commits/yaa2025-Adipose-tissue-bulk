# Adipose Bulk RNA‑seq Stress Analysis (DESeq2 + SVA)

This repo is a **GitHub‑ready** version of your analysis. It keeps your full pipeline (DESeq2, SVA, PCA, volcano, heatmap, GO) and saves outputs to `results/`.

> Groups: **HC = Control**, **C = CMVS**, **E = CSDS**.

---

## 🚀 Quick start

1) **Get the files** (this project):
   - Download the zip attached to our chat or clone your repo after you upload it.

2) **Install R packages** (first time only):
```r
source("requirements.R")
```

3) **Add your data**:
   - Put your raw count `.txt` files in `data/counts/`
   - Copy your metadata to `data/metadata.csv`
   - A template is provided at `data/metadata_template.csv`

4) **Run the pipeline**:
```r
source("analysis_pipeline.R")
# Example comparisons are at the bottom of the file (search for "USAGE").
```

Outputs are written to `results/`:
- `results/deseq2/` — CSVs
- `results/plots/`  — PCA/volcano/heatmap PDFs
- `results/go/`     — enrichment results

---

## 🧠 What’s inside (and why it’s longer than the “minimal” example)

Earlier I showed a **short, minimal script** just to illustrate the GitHub layout.  
This repo now includes your **full pipeline**, cleaned up for portability and reproducibility:
- No `setwd()`; all paths are relative to the project folder.
- A proper `requirements.R` that installs CRAN **and** Bioconductor packages.
- A `data/` folder with placeholders so you never commit big, private datasets by accident.
- `.gitignore` shields `results/` and big data files from Git.
- Functions kept: **SVA integration**, **run_deseq2_pipeline**, **PCA**, **volcano**, **DEG summary by time**, **heatmap**, **GO analysis**.

---

## 📁 Project tree

```
my-adipose-DESeq2/
│── README.md
│── LICENSE
│── .gitignore
│── requirements.R
│── analysis_pipeline.R
│── scripts/
│   └── run_examples.R
│── data/
│   ├── counts/
│   │   └── .gitkeep
│   ├── metadata_template.csv
│   └── README.md
│── results/
│   ├── .gitkeep
│   ├── deseq2/
│   ├── plots/
│   └── go/
```

---

## ⬆️ How to upload to GitHub (step‑by‑step)

### Option A — Web (easiest, no Terminal)
1. Go to https://github.com → **New repository** → name it `my-adipose-DESeq2` (Public/Private is your choice).
2. On your computer, **zip** this project folder (it’s already zipped in our chat download).
3. In your repo page, click **“Upload files”** → drag‑drop the entire zip contents (or the folder).
4. Add a commit message like “initial commit” → **Commit changes**.
5. Done! Your repo is live.

### Option B — Command line (git)
If you have git installed:

```bash
# 1) Go to the folder
cd path/to/my-adipose-DESeq2

# 2) Initialize and make the first commit
git init
git add .
git commit -m "Initial commit: RNA-seq DESeq2 + SVA pipeline"

# 3) Create an empty repo on GitHub (via the website). Copy its HTTPS URL, e.g.:
# https://github.com/<your-username>/my-adipose-DESeq2.git

# 4) Link and push
git branch -M main
git remote add origin https://github.com/<your-username>/my-adipose-DESeq2.git
git push -u origin main
```

> If you’re on Windows and want a GUI, **GitHub Desktop** works great. Choose “Add existing repository” and point it at this folder.

---

## 🧪 Data expectations

- **Counts**: one file per sample (`.txt`), 2 columns: `GeneID` and `Count`, **no header**.
- **Metadata**: CSV with at least these columns:
  - `Filename` (without `.txt`), `Group` ∈ {`Control`, `CMVS`, `CSDS`}, `Timepoint` ∈ {`zt_2`, `zt_6`, `zt_10`, `zt_14`, `zt_18`, `zt_22`}.

See `data/metadata_template.csv` and `data/README.md`.

---

## ❓ Need help

Ping me with any error message and I’ll point to the exact line and fix.
