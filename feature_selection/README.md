# Phase 1: Biological-Anchored Feature Selection Pipeline

This directory contains the reproducible workflow for condensing ~450,000 CpGs into the **198 core topological features** used by GT-Mamba. The pipeline follows a **Mechanism-Anchored Residual Ensemble (MARE)** strategy to ensure both biological interpretability and predictive robustness.

---

## 🛠️ Pipeline Overview

The selection process is organized into six sequential stages. Each script corresponds to a specific filtering or optimization task:

| Stage | Script Name | Core Task | Dimensionality |
| :--- | :--- | :--- | :--- |
| **01** | `01_partial_correlation_screening.py` | **Denoising**: Regressing out 6 major immune cell types to calculate pure age-related partial correlations. | 450,000 → 5,000 |
| **02** | `02_extract_candidate_matrix.py`      | **Alignment**: Merging 10 independent discovery cohorts into a unified candidate matrix for cross-cohort stability. | 5,000 (10 Cohorts) |
| **03** | `03_mechanism_anchoring_extraction.py`| **Anchoring**: Utilizing NCBI RefSeq coordinates to extract CpGs within regulatory windows of aging-related proteins. | Biological Priors |
| **04** | `04_mamba_ig_ranking.py`              | **Attribution**: Ranking sites via **Integrated Gradients (IG)** to quantify feature contribution in the Mamba architecture. | Feature Importance |
| **05** | `05_mamba_feature_sweep.py`           | **Scanning**: Iterative evaluation across feature scales (500-3000) to identify the optimal performance-complexity trade-off. | Scale Optimization |
| **06** | `06_stepwise_feature_refinement.py`   | **Locking**: Finalizing the 198-set via stepwise MAE optimization starting from the anchored baseline. | **Final 198 CpGs** |

---

## 🧬 Hierarchical Feature Tiers

The final **198 core CpGs** are structured into three hierarchical tiers to balance prior knowledge with data-driven discovery:

1. **Tier 1: Mandatory Anchors (121 CpGs)**
   CpG sites directly linked to the aging-related proteomic profile (ProtAge20). These sites act as biological "anchors" ensuring the model remains grounded in known aging mechanisms.

2. **Tier 2: Dual-Evidence Intersection (17 CpGs)**
   Sites validated by both biological prior knowledge (Genomic coordinates) and high partial correlation signals from the discovery datasets.

3. **Tier 3: Data-Driven Residuals (60 CpGs)**
   High-importance sites identified via AI attribution (Integrated Gradients). These sites capture essential aging signals that may lack explicit proteomic annotation but show strong topological consistency.

---

## 🧠 Why 198 CpGs?

Through our benchmarking, 198 sites provided the **optimal Pareto front** across three critical dimensions:

* **Biological Interpretability**: Over 70% of the signature (Tier 1 & 2) is explicitly anchored to aging-related pathways, ensuring the clock captures fundamental biological processes rather than statistical noise.
* **Computational Efficiency**: 198 is an optimal sequence length for the Graph Mamba block, ensuring the effective capture of genome-wide long-range dependencies while fully leveraging Mamba's hardware-aware linear complexity.
* **Cross-Platform Robustness**: These 198 sites show near-total overlap between **Illumina 450K** and **EPIC 850K** arrays, ensuring stable performance across different technologies.

---

## 🚀 Execution Guide

### Data Requirements
Ensure the `data/` directory contains the following reference files:
* `450k_annotation.csv` (Illumina Official Annotation)
* `ncbiRefSeq_hg19.csv` (Genomic Coordinates)
* `204aps.csv` / `20aps.csv` (Protein-gene mapping files)

### Running the Workflow
Execute the scripts in numerical order to replicate the full selection process:
```bash
python 01_partial_correlation_screening.py
python 02_extract_candidate_matrix.py
python 03_mechanism_anchoring_extraction.py
python 04_mamba_ig_ranking.py
python 05_mamba_feature_sweep.py
python 06_stepwise_feature_refinement.py
```