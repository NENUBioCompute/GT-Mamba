# Phase 2: Spatial Topology & Graph Modeling

This directory contains the logic and scripts for constructing the **Biologically-Informed Topology Graph** used by GT-Mamba. Unlike standard "all-to-all" attention mechanisms, GT-Mamba restricts information flow to functionally and spatially related CpG clusters to preserve biological inductive biases.

## 🕸️ Topology Construction Logic

The graph consists of **198 nodes** and is constructed through a multi-dimensional weighting strategy that integrates physical proximity and functional co-regulation:

### 1. Genomic Distance (Physical Proximity)
CpGs located within the same regulatory element (e.g., the same CpG island, enhancer, or promoter region) are more likely to exhibit coordinated methylation patterns. We utilize genomic coordinates (GRCh37/hg19) to link neighboring sites.

### 2. Co-methylation Correlation (Functional Linkage)
Using the large-scale discovery cohorts, we calculate the Pearson correlation coefficients to capture functional co-regulation during the aging process.

### 3. K-Nearest Neighbor (KNN) Sparsification
To maintain high computational efficiency, a **K=5** sparsification strategy is applied. Each CpG node is connected only to its top 5 most significant "neighbors," ranked by a fusion score of physical distance and functional correlation.

---

## 📁 Key File Outputs

Executing the construction script will generate the following artifacts in the `data/graph/` directory:

* **`methylation_graph.pt`**: A PyTorch Geometric (PyG) data object containing the `edge_index`. This file is loaded by GT-Mamba during initialization to guide the masked attention mechanism.
* **`methylation_graph.graphml`**: An XML-based graph format designed for visualization in software such as **Gephi** or **Cytoscape**.
* **`graph_summary.csv`**: A statistical summary detailing the structural properties (e.g., degree, chromosome location) of the 198-CpG nodes.

---

## 🛠️ Usage

To construct the topology graph, simply execute the pipeline script:

```bash
# Standard construction (K=5)
python build_topology_graph.py

# Custom construction (e.g., experimenting with K=10)
python build_topology_graph.py --k 10
```

## 🎨 Visualization
The topological representation of the 198-CpG clock reveals distinct functional clusters centered around master aging-regulatory genes such as *ELOVL2*, *FHL2*, and *KLF4*.