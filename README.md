
# GT-Mamba: A Topology-Aware Graph-State Space Model for Robust and Interpretable Epigenetic Age Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Framework](https://img.shields.io/badge/PyTorch-2.0%2B-red)](https://pytorch.org/)

This repository contains the official implementation of the paper:  
**GT-Mamba: A Topology-Aware Graph-State Space Model for Robust and Interpretable Epigenetic Age Prediction**

> **Abstract:** We propose **GT-Mamba**, a novel architecture that integrates a **Structure-Aware Graph Transformer** with the **Mamba** state space model. By explicitly modeling CpG topological correlations and genome-wide long-range dependencies, GT-Mamba achieves state-of-the-art precision and superior robustness across heterogeneous cohorts.

![Model Architecture](Figure1.png)
*(Figure 1: Overview of the GT-Mamba Architecture.)*

## 📂 Directory Structure

The codebase is organized as follows:

```text
GT-Mamba/
├── checkpoints/          # Pre-trained model weights
│   └── best_model(5).pth # ✅ Best performing model checkpoint
├── data/                 # Data directory
│   └── methylation_matrix_198.csv  # ✅ Core dataset (198 CpGs)
├── data_preparation/     # Scripts for data preprocessing
│   ├── impute_methylation.R
│   └── run_preprocess.py
├── feature_process/      # Feature selection and processing modules
├── graph_construction/   # Graph topology generation
│   └── build-graph.py
├── model_training/       # Main model training scripts
│   └── train-graph.py
├── requirements.txt      # Python dependencies
└── README.md             # Documentation

```

## 🛠️ Installation

We recommend using **Anaconda** to manage the environment, as `mamba-ssm` requires specific CUDA setups.

### 1. Create Environment

```bash
conda create -n gtmamba python=3.10
conda activate gtmamba

```

### 2. Install PyTorch

Please install PyTorch matching your CUDA version (e.g., CUDA 11.8):

```bash
pip install torch==2.0.1 --index-url [https://download.pytorch.org/whl/cu118](https://download.pytorch.org/whl/cu118)

```

### 3. Install Mamba & Graph Dependencies

```bash
# Install Graph Neural Network libraries
pip install torch_geometric

# Install Mamba (mambapy) and other dependencies
# Note: Ensure mamba-ssm or mambapy is installed correctly for your CUDA version
pip install -r requirements.txt

```

## 🧬 Data Availability

* **Core Dataset (Included):** The processed dataset containing the 198 core CpG sites is provided in `data/methylation_matrix_198.csv`.
* **Public Datasets:** Raw external validation cohorts (e.g., GSE40279, GSE61496) are publicly available on NCBI GEO.

## 🚀 Usage

### Step 1: Preprocessing

You can use the provided Python wrapper or run the R script directly for imputation:

```bash
# Option A: Python wrapper
python data_preparation/run_preprocess.py

# Option B: Direct R script
Rscript data_preparation/impute_methylation.R

```

### Step 2: Graph Construction

Construct the graph topology based on the data:

```bash
python graph_construction/build-graph.py --data data/methylation_matrix_198.csv --k 5

```

### Step 3: Model Training

Train the GT-Mamba model using the graph data:

```bash
python model_training/train-graph.py --mode train \
                                     --data_path data/methylation_matrix_198.csv \
                                     --epochs 400 \
                                     --batch_size 8 \
                                     --save_dir checkpoints/

```

### Prediction

To predict biological age using the provided pre-trained weights (`best_model(5).pth`):

```bash
python model_training/train-graph.py --mode predict \
                                     --data_path data/new_samples.csv \
                                     --checkpoint checkpoints/best_model(5).pth \
                                     --output results.csv

```

## 🔗 Citation

If you use this code, please cite our paper:

```bibtex
@article{Wang2026GTMamba,
  title={GT-Mamba: A Topology-Aware Graph-State Space Model for Robust and Interpretable Epigenetic Age Prediction},
  author={Wang, Han and Wang, Hui and Tong, Yanting and Liu, Yuanyuan and Jing, Qu and Zhang, Li},
  journal={Journal Title Here},
  year={2026},
  note={Under Review}
}

```

## 📧 Contact

* **Li Zhang**: [lizhang@ccut.edu.cn](mailto:lizhang@ccut.edu.cn)
* **Han Wang**: [wanghan101@nenu.edu.cn](mailto:wanghan101@nenu.edu.cn)

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

```

```
