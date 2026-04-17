# graph_modeling/build_topology_graph.py
import os
import argparse
import warnings
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import torch
from torch_geometric.data import Data
import networkx as nx
from tqdm import tqdm

warnings.filterwarnings('ignore')


# ==============================================================================
# GT-Mamba Graph Construction Module (Phase 2)
# Purpose: Construct a biologically meaningful topology graph for the 198 core CpGs.
# Mechanism: Edges are weighted based on a geometric integration of:
#            1. Co-methylation Correlation (Pearson)
#            2. Genomic Distance (kb/Mb decay)
#            3. Regulatory Regions (CpG Island/Shore/Shelf weighting)
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="GT-Mamba Topology Graph Builder")
    parser.add_argument("--k", type=int, default=5, help="KNN sparsification parameter (default: 5)")
    parser.add_argument("--alpha", type=float, default=0.55, help="Baseline regulatory weight (default: 0.55)")
    return parser.parse_args()


def main():
    args = parse_args()
    K_VAL = args.k
    ALPHA_VAL = args.alpha

    print(f"🚀 Starting Graph Construction (K={K_VAL}, Alpha={ALPHA_VAL})")

    # 1. Automatically locate the project root and define paths
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = os.path.dirname(CURRENT_DIR)

    # Input file paths
    DATA_DIR = os.path.join(BASE_DIR, "data")
    ANNO_FILE = os.path.join(DATA_DIR, "450k_annotation.csv")
    REFSEQ_FILE = os.path.join(DATA_DIR, "ncbiRefSeq_hg19.csv")
    MATRIX_FILE = os.path.join(DATA_DIR, "methylation_matrix_198.csv")

    # Output file paths
    GRAPH_OUT_DIR = os.path.join(DATA_DIR, "graph")
    os.makedirs(GRAPH_OUT_DIR, exist_ok=True)

    print("📥 Loading data...")
    annotation_df = pd.read_csv(ANNO_FILE)
    refseq_df = pd.read_csv(REFSEQ_FILE)
    df = pd.read_csv(MATRIX_FILE, index_col=0)

    # Clean columns and extract methylation data matrix
    df.columns = df.columns.str.strip()
    methylation_data = df.drop(columns=['Age', 'Split', 'ID', 'Dataset'], errors='ignore')
    n_cpgs = len(methylation_data.columns)
    n_samples = len(methylation_data)
    print(f"✅ Methylation matrix loaded: {n_samples} samples, {n_cpgs} CpGs")

    # ==================== 2. Integrate Biological Annotations ====================
    print("\n🧬 Integrating biological annotations...")
    cpg_id_column = 'Name'
    annotation_df[cpg_id_column] = annotation_df[cpg_id_column].astype(str).str.strip()
    available_cpgs = [cpg.strip() for cpg in methylation_data.columns]

    # Extract relevant positions and regulatory regions
    cpg_positions = annotation_df[annotation_df[cpg_id_column].isin(available_cpgs)].copy()
    cpg_positions = cpg_positions[['Name', 'chr', 'pos', 'Relation_to_Island', 'UCSC_RefGene_Name']].copy()
    cpg_positions = cpg_positions.rename(columns={
        'Name': 'CpG_ID', 'chr': 'CpG_Chromosome', 'pos': 'Position',
        'Relation_to_Island': 'Region', 'UCSC_RefGene_Name': 'Gene'
    })

    # Standardize RefSeq data
    refseq_df.columns = refseq_df.columns.str.strip()
    gene_regions = refseq_df[['name2', 'chrom', 'txStart', 'txEnd']].copy().rename(columns={
        'name2': 'Gene', 'chrom': 'Gene_Chromosome', 'txStart': 'Start', 'txEnd': 'End'
    })
    merged_df = pd.merge(cpg_positions, gene_regions, on='Gene', how='left')

    # ==================== 3. CpG Node Class Definition ====================
    class CpGSiteInfo:
        def __init__(self, cpg_id):
            self.cpg_id = cpg_id
            self.has_position = False
            self.has_regulatory = False
            self.chromosome = None
            self.position = None
            self.regulatory_weight = ALPHA_VAL  # Dynamic baseline weight
            self.site_type = None

        def set_position_info(self, chromosome, position):
            if pd.notna(chromosome) and pd.notna(position):
                self.chromosome = str(chromosome).strip()
                if not self.chromosome.startswith('chr'):
                    self.chromosome = f'chr{self.chromosome}'
                self.position = int(position)
                self.has_position = True

        def set_regulatory_info(self, region, gene_info):
            if pd.notna(region):
                region_str = str(region).upper()
                if any(k in region_str for k in ['ISLAND', 'TSST1500']):
                    self.regulatory_weight = 1.0
                    self.has_regulatory = True
                elif any(k in region_str for k in ['SHORE', 'SHELF']):
                    self.regulatory_weight = 0.7
                    self.has_regulatory = True
                elif 'OPENSEA' in region_str:
                    self.regulatory_weight = 0.4
                    self.has_regulatory = True

            # Give a slight boost if associated with a known gene
            if pd.notna(gene_info) and self.has_regulatory:
                self.regulatory_weight = min(1.0, self.regulatory_weight + 0.1)

        def determine_site_type(self):
            if self.has_position and self.has_regulatory:
                self.site_type = 'A'
            elif not self.has_position and self.has_regulatory:
                self.site_type = 'B'
            elif self.has_position and not self.has_regulatory:
                self.site_type = 'C'
            else:
                self.site_type = 'D'
            return self.site_type

    # Initialize all nodes
    cpg_info_dict = {}
    for cpg in methylation_data.columns:
        info = CpGSiteInfo(cpg)
        cpg_data = merged_df[merged_df['CpG_ID'] == cpg]
        if not cpg_data.empty:
            info.set_position_info(cpg_data.iloc[0]['CpG_Chromosome'], cpg_data.iloc[0]['Position'])
            info.set_regulatory_info(cpg_data.iloc[0]['Region'], cpg_data.iloc[0]['Gene'])
        info.determine_site_type()
        cpg_info_dict[cpg] = info

    # ==================== 4. Multi-modal Adjacency Weights Calculation ====================
    print("\n🧮 Calculating Multi-modal Adjacency Weights...")

    def calculate_correlation_matrix(data_matrix):
        """Calculate Pearson correlation for functional linkage."""
        n_features = data_matrix.shape[1]
        corr_matrix = np.ones((n_features, n_features))
        data_array = data_matrix.values.T
        for i in tqdm(range(n_features), desc="Pearson Correlation"):
            for j in range(i + 1, n_features):
                x, y = data_array[i], data_array[j]
                mask = ~(np.isnan(x) | np.isnan(y))
                if mask.sum() < 3:
                    corr = 0.0
                else:
                    corr, _ = pearsonr(x[mask], y[mask])
                    corr = 0.0 if np.isnan(corr) else abs(corr)
                corr_matrix[i, j] = corr_matrix[j, i] = corr
        return corr_matrix

    def calculate_distance_weight(distance):
        """Calculate weight decay based on genomic distance."""
        if distance < 1000:
            return np.exp(-distance / 500)
        elif distance < 10000:
            return 0.9 * np.exp(-(distance - 1000) / 5000)
        elif distance < 100000:
            return 0.6 * np.exp(-(distance - 10000) / 30000)
        elif distance < 1000000:
            return 0.3 * np.exp(-(distance - 100000) / 300000)
        else:
            return 0.1 * np.exp(-(distance - 1000000) / 1000000)

    def build_dist_and_reg_matrices(cpg_list):
        """Construct distance and regulatory matrices."""
        n = len(cpg_list)
        dist_matrix = np.ones((n, n))
        reg_matrix = np.ones((n, n))

        # Pre-calculate regulatory geometry (geometric mean of node weights)
        reg_weights = np.array([cpg_info_dict[c].regulatory_weight for c in cpg_list])
        for i in range(n):
            for j in range(n):
                reg_matrix[i, j] = np.sqrt(reg_weights[i] * reg_weights[j])

        # Distance calculation
        for i in tqdm(range(n), desc="Genomic Distance"):
            for j in range(i + 1, n):
                info_i, info_j = cpg_info_dict[cpg_list[i]], cpg_info_dict[cpg_list[j]]
                if info_i.has_position and info_j.has_position:
                    if info_i.chromosome == info_j.chromosome:
                        weight = calculate_distance_weight(abs(info_i.position - info_j.position))
                    else:
                        weight = 0.001  # Severe penalty for inter-chromosomal links
                else:
                    weight = 1.0  # Neutral if location is unknown
                dist_matrix[i, j] = dist_matrix[j, i] = weight

        return dist_matrix, reg_matrix

    cpg_list = list(methylation_data.columns)
    corr_matrix = calculate_correlation_matrix(methylation_data)
    dist_matrix, reg_matrix = build_dist_and_reg_matrices(cpg_list)

    # Combine Weights (Geometric Mean Integration)
    print("\n🔗 Fusing Modalities into Final Graph...")
    combined_matrix = np.zeros((n_cpgs, n_cpgs))
    for i in range(n_cpgs):
        for j in range(i + 1, n_cpgs):
            info_i, info_j = cpg_info_dict[cpg_list[i]], cpg_info_dict[cpg_list[j]]
            factors = [corr_matrix[i, j]]

            # Incorporate Distance Penalty
            if info_i.has_position and info_j.has_position:
                factors.append(dist_matrix[i, j])
            else:
                factors.append(1.0)

            # Incorporate Regulatory Knowledge
            if info_i.has_regulatory or info_j.has_regulatory:
                factors.append(reg_matrix[i, j])
            else:
                factors.append(ALPHA_VAL)

            # Compute the geometric mean of all multi-modal factors
            w = np.prod(factors) ** (1 / len(factors))
            combined_matrix[i, j] = combined_matrix[j, i] = w

    # ==================== 5. KNN Sparsification ====================
    print(f"\n✂️ Sparsifying Graph (KNN-{K_VAL})...")
    sparse_adj = np.zeros_like(combined_matrix)
    for i in range(n_cpgs):
        w = combined_matrix[i].copy()
        w[i] = -np.inf  # Exclude self-loops during KNN selection
        top_k = np.argsort(w)[-K_VAL:]
        sparse_adj[i, top_k] = combined_matrix[i, top_k]

    # Ensure graph symmetry for undirected graph operations
    sparse_adj = np.maximum(sparse_adj, sparse_adj.T)

    # ==================== 6. Export PyG and NetworkX Artifacts ====================
    print("\n💾 Exporting Graph Artifacts (Standardized Filenames)...")
    edges = np.nonzero(sparse_adj)
    edge_index = torch.tensor(edges, dtype=torch.long)
    node_features = torch.tensor(methylation_data.values.T, dtype=torch.float)
    edge_weights = torch.tensor(sparse_adj[edges], dtype=torch.float).unsqueeze(1)

    # Standardized output files without suffix for streamlined downstream loading
    pt_file = os.path.join(GRAPH_OUT_DIR, 'methylation_graph.pt')
    graphml_file = os.path.join(GRAPH_OUT_DIR, 'methylation_graph.graphml')
    csv_file = os.path.join(GRAPH_OUT_DIR, 'graph_summary.csv')

    # Save PyTorch Geometric Data object
    graph_data = Data(x=node_features, edge_index=edge_index, edge_attr=edge_weights)
    torch.save(graph_data, pt_file)
    print(
        f"  ✓ PyG Graph saved -> {os.path.basename(pt_file)} (Nodes: {graph_data.num_nodes}, Edges: {graph_data.num_edges})")

    # Construct NetworkX Graph for visualization and structural analysis
    G = nx.Graph()
    for i, cpg in enumerate(cpg_list):
        G.add_node(i, cpg_id=str(cpg))
    for i in range(edges[0].shape[0]):
        u, v = edges[0][i], edges[1][i]
        if u != v:
            G.add_edge(u, v, weight=float(sparse_adj[u, v]))
    nx.write_graphml(G, graphml_file)

    # Save Graph structural summary
    summary_df = pd.DataFrame([{
        'node_id': i, 'cpg_id': cpg,
        'chromosome': cpg_info_dict[cpg].chromosome if cpg_info_dict[cpg].has_position else "Unknown",
        'degree': len(list(G.neighbors(i))),
        'clustering_coefficient': nx.clustering(G, weight='weight').get(i, 0.0)
    } for i, cpg in enumerate(cpg_list)])

    summary_df.to_csv(csv_file, index=False)
    print(f"  ✓ Graph Summary CSV saved -> {os.path.basename(csv_file)}")
    print("\n🎉 Topology Construction Complete!")


if __name__ == "__main__":
    main()