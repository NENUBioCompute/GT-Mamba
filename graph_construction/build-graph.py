# import pandas as pd
# import numpy as np
# from scipy.stats import pearsonr
# import torch
# from torch_geometric.data import Data
# import networkx as nx
# from tqdm import tqdm
# import warnings
#
# warnings.filterwarnings('ignore')
#
# # ==================== 1. 加载数据 ====================
# print("正在加载数据...")
# annotation_df = pd.read_csv("/datapool/home/info_wang/wanghui/file/450k_annotation.csv")
# refseq_df = pd.read_csv("/datapool/home/info_wang/wanghui/file/ncbiRefSeq_hg19.csv")
# df = pd.read_csv("/datapool/home/info_wang/wanghui/file/CPG/methylation_matrix_198.csv", index_col=0)
#
# # 清理列名
# df.columns = df.columns.str.strip()
# methylation_data = df.drop(columns=['Age', 'Split', 'ID'])
# n_cpgs = len(methylation_data.columns)
# n_samples = len(methylation_data)
#
# print(f"甲基化数据: {n_samples}个样本, {n_cpgs}个CpG位点")
#
# # ==================== 2. 整合注释信息 ====================
# print("\n=== 整合注释信息 ===")
# cpg_id_column = 'Name'
#
# # 清理位点名称
# annotation_df[cpg_id_column] = annotation_df[cpg_id_column].astype(str).str.strip()
# available_cpgs = [cpg.strip() for cpg in methylation_data.columns]
#
# # 筛选匹配的位点
# cpg_positions = annotation_df[annotation_df[cpg_id_column].isin(available_cpgs)].copy()
# print(f"在注释文件中找到 {len(cpg_positions)} 个匹配的CpG位点")
#
# # 选择需要的列并重命名
# cpg_positions = cpg_positions[['Name', 'chr', 'pos', 'Relation_to_Island', 'UCSC_RefGene_Name']].copy()
# cpg_positions = cpg_positions.rename(columns={
#     'Name': 'CpG_ID',
#     'chr': 'CpG_Chromosome',
#     'pos': 'Position',
#     'Relation_to_Island': 'Region',
#     'UCSC_RefGene_Name': 'Gene'
# })
#
# # ==================== 3. 处理基因区域信息 ====================
# refseq_df.columns = refseq_df.columns.str.strip()
# gene_regions = refseq_df[['name2', 'chrom', 'txStart', 'txEnd']].copy()
# gene_regions = gene_regions.rename(columns={
#     'name2': 'Gene',
#     'chrom': 'Gene_Chromosome',
#     'txStart': 'Start',
#     'txEnd': 'End'
# })
#
# # 合并注释信息和基因区域信息
# merged_df = pd.merge(cpg_positions, gene_regions, on='Gene', how='left')
#
#
# # ==================== 4. 创建CpG位点信息类 ====================
# class CpGSiteInfo:
#     """存储每个CpG位点的完整信息"""
#
#     def __init__(self, cpg_id):
#         self.cpg_id = cpg_id
#         self.has_position = False
#         self.has_regulatory = False
#         self.chromosome = None
#         self.position = None
#         self.regulatory_weight = None
#         self.site_type = None  # 'A', 'B', 'C', 'D'
#
#     def set_position_info(self, chromosome, position):
#         """设置位置信息"""
#         if pd.notna(chromosome) and pd.notna(position):
#             self.chromosome = str(chromosome).strip()
#             if not self.chromosome.startswith('chr'):
#                 self.chromosome = f'chr{self.chromosome}'
#             self.position = int(position)
#             self.has_position = True
#
#     def set_regulatory_info(self, region, gene_info):
#         """根据注释信息设置调控权重"""
#         # 默认权重
#         self.regulatory_weight = 0.55  # 保守中间值
#         self.has_regulatory = False
#
#         # 如果有区域信息
#         if pd.notna(region):
#             region_str = str(region).upper()
#
#             # 第1级：核心调控区域 (权重=1.0)
#             if any(keyword in region_str for keyword in ['ISLAND', 'TSST1500']):
#                 self.regulatory_weight = 1.0
#                 self.has_regulatory = True
#
#             # 第2级：潜在调控区域 (权重=0.7)
#             elif any(keyword in region_str for keyword in ['SHORE', 'SHELF']):
#                 self.regulatory_weight = 0.7
#                 self.has_regulatory = True
#
#             # 第3级：非调控区域 (权重=0.4)
#             elif 'OPENSEA' in region_str:
#                 self.regulatory_weight = 0.4
#                 self.has_regulatory = True
#
#         # 如果有基因信息，可以提升权重
#         if pd.notna(gene_info) and self.has_regulatory:
#             # 如果有基因关联，稍微提升权重
#             self.regulatory_weight = min(1.0, self.regulatory_weight + 0.1)
#
#     def determine_site_type(self):
#         """确定位点类型"""
#         if self.has_position and self.has_regulatory:
#             self.site_type = 'A'
#         elif not self.has_position and self.has_regulatory:
#             self.site_type = 'B'
#         elif self.has_position and not self.has_regulatory:
#             self.site_type = 'C'
#         else:  # not has_position and not has_regulatory
#             self.site_type = 'D'
#         return self.site_type
#
#
# # ==================== 5. 初始化所有CpG位点信息 ====================
# print("\n=== 初始化CpG位点信息 ===")
# cpg_info_dict = {}
#
# for cpg in methylation_data.columns:
#     info = CpGSiteInfo(cpg)
#
#     # 查找位置信息
#     cpg_data = merged_df[merged_df['CpG_ID'] == cpg]
#     if not cpg_data.empty:
#         info.set_position_info(
#             cpg_data.iloc[0]['CpG_Chromosome'],
#             cpg_data.iloc[0]['Position']
#         )
#
#         # 设置调控信息
#         info.set_regulatory_info(
#             cpg_data.iloc[0]['Region'],
#             cpg_data.iloc[0]['Gene']
#         )
#
#     # 确定位点类型
#     info.determine_site_type()
#     cpg_info_dict[cpg] = info
#
# # 统计位点类型
# type_counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
# for info in cpg_info_dict.values():
#     type_counts[info.site_type] += 1
#
# print(f"位点类型统计:")
# print(f"  A类(有位置+有调控): {type_counts['A']} ({type_counts['A'] / n_cpgs * 100:.1f}%)")
# print(f"  B类(无位置+有调控): {type_counts['B']} ({type_counts['B'] / n_cpgs * 100:.1f}%)")
# print(f"  C类(有位置+无调控): {type_counts['C']} ({type_counts['C'] / n_cpgs * 100:.1f}%)")
# print(f"  D类(无位置+无调控): {type_counts['D']} ({type_counts['D'] / n_cpgs * 100:.1f}%)")
#
# # ==================== 6. 计算相关性权重矩阵 ====================
# print("\n=== 计算相关性权重矩阵 ===")
#
#
# def calculate_correlation_matrix_optimized(data_matrix):
#     """优化版本的相关性计算"""
#     n_features = data_matrix.shape[1]
#     corr_matrix = np.ones((n_features, n_features))
#
#     # 转换为numpy数组并转置为(特征数, 样本数)
#     data_array = data_matrix.values.T
#
#     print("计算皮尔森相关系数...")
#     for i in tqdm(range(n_features), desc="处理位点"):
#         for j in range(i + 1, n_features):
#             # 获取两个位点的数据
#             x = data_array[i]
#             y = data_array[j]
#
#             # 去除缺失值
#             mask = ~(np.isnan(x) | np.isnan(y))
#             x_valid = x[mask]
#             y_valid = y[mask]
#
#             # 至少需要3个有效样本
#             if len(x_valid) < 3:
#                 corr = 0.0
#             else:
#                 # 计算皮尔森相关系数
#                 corr, _ = pearsonr(x_valid, y_valid)
#                 if np.isnan(corr):
#                     corr = 0.0
#
#             # 存储到矩阵（使用绝对值）
#             corr_matrix[i, j] = abs(corr)
#             corr_matrix[j, i] = abs(corr)
#
#     return corr_matrix
#
#
# # 计算相关性矩阵
# corr_matrix = calculate_correlation_matrix_optimized(methylation_data)
# print(f"相关性矩阵计算完成，形状: {corr_matrix.shape}")
#
# # ==================== 7. 计算距离权重矩阵 ====================
# print("\n=== 计算距离权重矩阵 ===")
#
#
# def calculate_distance_weight(distance):
#     """分段计算距离权重"""
#     if distance < 1000:  # < 1kb
#         return np.exp(-distance / 500)
#     elif distance < 10000:  # 1kb - 10kb
#         return 0.9 * np.exp(-(distance - 1000) / 5000)
#     elif distance < 100000:  # 10kb - 100kb
#         return 0.6 * np.exp(-(distance - 10000) / 30000)
#     elif distance < 1000000:  # 100kb - 1Mb
#         return 0.3 * np.exp(-(distance - 100000) / 300000)
#     else:  # > 1Mb
#         return 0.1 * np.exp(-(distance - 1000000) / 1000000)
#
#
# def calculate_distance_weight_matrix(cpg_info_dict, cpg_list):
#     """计算距离权重矩阵"""
#     n = len(cpg_list)
#     dist_weight_matrix = np.ones((n, n))  # 默认值为1（中性）
#
#     # 创建位点索引映射
#     cpg_to_idx = {cpg: i for i, cpg in enumerate(cpg_list)}
#
#     # 只对有位置信息的位点进行计算
#     valid_cpgs = [cpg for cpg in cpg_list if cpg_info_dict[cpg].has_position]
#     print(f"有位置信息的位点数量: {len(valid_cpgs)}/{n}")
#
#     # 计算同染色体位点的距离
#     chromosome_groups = {}
#     for cpg in valid_cpgs:
#         info = cpg_info_dict[cpg]
#         chrom = info.chromosome
#         if chrom not in chromosome_groups:
#             chromosome_groups[chrom] = []
#         chromosome_groups[chrom].append((cpg, info.position, cpg_to_idx[cpg]))
#
#     print(f"分组到 {len(chromosome_groups)} 条染色体")
#
#     # 对每条染色体内部的位点计算距离
#     for chrom, sites in tqdm(chromosome_groups.items(), desc="计算染色体内部距离"):
#         if len(sites) < 2:
#             continue
#
#         # 按位置排序
#         sites.sort(key=lambda x: x[1])
#
#         for i in range(len(sites)):
#             cpg1, pos1, idx1 = sites[i]
#             for j in range(i + 1, len(sites)):
#                 cpg2, pos2, idx2 = sites[j]
#                 distance = abs(pos1 - pos2)
#                 weight = calculate_distance_weight(distance)
#                 dist_weight_matrix[idx1, idx2] = weight
#                 dist_weight_matrix[idx2, idx1] = weight
#
#     # 不同染色体的位点对，权重设为0.001
#     different_chromosome_pairs = 0
#     for i in range(n):
#         for j in range(i + 1, n):
#             info_i = cpg_info_dict[cpg_list[i]]
#             info_j = cpg_info_dict[cpg_list[j]]
#
#             if info_i.has_position and info_j.has_position:
#                 if info_i.chromosome != info_j.chromosome:
#                     dist_weight_matrix[i, j] = 0.001
#                     dist_weight_matrix[j, i] = 0.001
#                     different_chromosome_pairs += 1
#
#     print(f"不同染色体位点对数量: {different_chromosome_pairs}")
#     return dist_weight_matrix
#
#
# # 计算距离权重矩阵
# cpg_list = list(methylation_data.columns)
# dist_weight_matrix = calculate_distance_weight_matrix(cpg_info_dict, cpg_list)
#
# # ==================== 8. 计算调控权重矩阵 ====================
# print("\n=== 计算调控权重矩阵 ===")
#
#
# def calculate_regulatory_weight_matrix(cpg_info_dict, cpg_list):
#     """计算调控权重矩阵"""
#     n = len(cpg_list)
#     reg_weight_matrix = np.ones((n, n))
#
#     # 获取每个位点的调控权重
#     reg_weights = np.array([cpg_info_dict[cpg].regulatory_weight for cpg in cpg_list])
#
#     # 计算调控权重矩阵：使用几何平均 R_ij = sqrt(r_i * r_j)
#     for i in range(n):
#         for j in range(n):
#             if i == j:
#                 reg_weight_matrix[i, j] = reg_weights[i]
#             else:
#                 reg_weight_matrix[i, j] = np.sqrt(reg_weights[i] * reg_weights[j])
#
#     return reg_weight_matrix
#
#
# # 计算调控权重矩阵
# reg_weight_matrix = calculate_regulatory_weight_matrix(cpg_info_dict, cpg_list)
#
# # ==================== 9. 计算综合权重矩阵 ====================
# print("\n=== 计算综合权重矩阵 ===")
#
#
# def calculate_combined_weight(corr_weight, dist_weight, reg_weight,
#                               info1, info2):
#     """根据信息完整性动态计算综合权重"""
#
#     # 确定有效因子数
#     n_factors = 1  # 相关性总是有
#     factors = [corr_weight]
#
#     # 距离因子
#     if info1.has_position and info2.has_position:
#         n_factors += 1
#         factors.append(dist_weight)
#     else:
#         # 无位置信息时，距离权重设为1（中性）
#         factors.append(1.0)
#
#     # 调控因子
#     if info1.has_regulatory or info2.has_regulatory:
#         n_factors += 1
#         factors.append(reg_weight)
#     else:
#         # 两者都无调控信息时，调控权重设为0.55（保守值）
#         factors.append(0.55)
#
#     # 计算几何平均
#     if n_factors == 3:
#         return np.prod(factors) ** (1 / 3)
#     elif n_factors == 2:
#         # 只使用相关性和调控（或无位置但有调控的情况）
#         return (corr_weight * factors[2]) ** (1 / 2)
#     else:  # n_factors == 1
#         return corr_weight
#
#
# def calculate_combined_weight_matrix(corr_matrix, dist_matrix, reg_matrix,
#                                      cpg_info_dict, cpg_list):
#     """计算综合权重矩阵"""
#     n = len(cpg_list)
#     combined_matrix = np.zeros((n, n))
#     info_completeness = np.zeros((n, n), dtype='U3')  # 存储信息完整性
#
#     print("计算综合权重...")
#     for i in tqdm(range(n), desc="处理位点"):
#         info_i = cpg_info_dict[cpg_list[i]]
#         for j in range(i, n):
#             info_j = cpg_info_dict[cpg_list[j]]
#
#             # 计算综合权重
#             weight = calculate_combined_weight(
#                 corr_matrix[i, j],
#                 dist_matrix[i, j],
#                 reg_matrix[i, j],
#                 info_i,
#                 info_j
#             )
#
#             combined_matrix[i, j] = weight
#             combined_matrix[j, i] = weight
#
#             # 记录信息完整性
#             completeness = f"{info_i.site_type}{info_j.site_type}"
#             info_completeness[i, j] = completeness
#             info_completeness[j, i] = completeness
#
#     return combined_matrix, info_completeness
#
#
# # 计算综合权重矩阵
# combined_matrix, info_completeness = calculate_combined_weight_matrix(
#     corr_matrix, dist_weight_matrix, reg_weight_matrix,
#     cpg_info_dict, cpg_list
# )
#
# # ==================== 10. 稀疏化处理 ====================
# print("\n=== 稀疏化处理 ===")
#
#
# def sparse_top_k(matrix, k=5):
#     """每个节点只保留权重最高的k个连接"""
#     n = matrix.shape[0]
#     sparse_matrix = np.zeros_like(matrix)
#
#     for i in range(n):
#         # 获取当前节点的所有权重
#         weights = matrix[i].copy()
#         weights[i] = -np.inf  # 排除自身
#
#         # 找出top-k个邻居
#         top_k_indices = np.argsort(weights)[-k:]
#         sparse_matrix[i, top_k_indices] = matrix[i, top_k_indices]
#
#     # 确保对称性
#     sparse_matrix = np.maximum(sparse_matrix, sparse_matrix.T)
#
#     # 计算稀疏化后的统计
#     edge_count = np.sum(sparse_matrix > 0) // 2
#     print(f"稀疏化后边数量: {edge_count}")
#     print(f"平均每个节点连接数: {2 * edge_count / n:.2f}")
#
#     return sparse_matrix
#
#
# # 稀疏化处理
# k_neighbors = 5
# sparse_adj_matrix = sparse_top_k(combined_matrix, k=k_neighbors)
#
# # ==================== 11. 构建图结构 ====================
# print("\n=== 构建图结构 ===")
#
# # 提取边
# edges = np.nonzero(sparse_adj_matrix)
# edge_index = torch.tensor(edges, dtype=torch.long)
#
# # 节点特征：甲基化数据转置
# node_features = torch.tensor(methylation_data.values.T, dtype=torch.float)
#
# # 边权重
# edge_weights = torch.tensor(sparse_adj_matrix[edges], dtype=torch.float)
#
# # 创建PyG图数据
# graph_data = Data(
#     x=node_features,
#     edge_index=edge_index,
#     edge_attr=edge_weights.unsqueeze(1)  # 转换为列向量
# )
#
# print(f"图数据信息:")
# print(f"- 节点数: {graph_data.num_nodes}")
# print(f"- 边数: {graph_data.num_edges}")
# print(f"- 节点特征维度: {graph_data.num_node_features}")
# print(f"- 边特征维度: {graph_data.num_edge_features}")
#
# # ==================== 12. 创建NetworkX图（用于分析） ====================
# print("\n=== 创建NetworkX图 ===")
# G = nx.Graph()
#
# # 添加节点
# for i, cpg in enumerate(cpg_list):
#     info = cpg_info_dict[cpg]
#     G.add_node(i,
#                cpg_id=cpg,
#                site_type=info.site_type,
#                chromosome=info.chromosome if info.has_position else None,
#                position=info.position if info.has_position else None,
#                regulatory_weight=info.regulatory_weight)
#
# # 添加边
# for i in range(edges[0].shape[0]):
#     source = edges[0][i]
#     target = edges[1][i]
#     if source != target:
#         weight = sparse_adj_matrix[source, target]
#         # 获取信息完整性
#         completeness = info_completeness[source, target]
#         G.add_edge(source, target,
#                    weight=weight,
#                    info_completeness=completeness,
#                    correlation=corr_matrix[source, target],
#                    distance_weight=dist_weight_matrix[source, target],
#                    regulatory_weight=reg_weight_matrix[source, target])
#
# print(f"NetworkX图信息:")
# print(f"- 节点数: {G.number_of_nodes()}")
# print(f"- 边数: {G.number_of_edges()}")
# print(f"- 图密度: {nx.density(G):.6f}")
#
# # ==================== 13. 图统计分析 ====================
# print("\n=== 图统计分析 ===")
#
# # 按位点类型统计连接
# type_connections = {}
# for u, v, data in G.edges(data=True):
#     u_type = G.nodes[u]['site_type']
#     v_type = G.nodes[v]['site_type']
#     edge_type = f"{u_type}-{v_type}"
#
#     if edge_type not in type_connections:
#         type_connections[edge_type] = {
#             'count': 0,
#             'total_weight': 0,
#             'avg_correlation': 0,
#             'avg_distance_weight': 0
#         }
#
#     type_connections[edge_type]['count'] += 1
#     type_connections[edge_type]['total_weight'] += data['weight']
#     type_connections[edge_type]['avg_correlation'] += data['correlation']
#     type_connections[edge_type]['avg_distance_weight'] += data['distance_weight']
#
# print("按连接类型统计:")
# for edge_type, stats in sorted(type_connections.items()):
#     count = stats['count']
#     avg_weight = stats['total_weight'] / count
#     avg_corr = stats['avg_correlation'] / count
#     avg_dist = stats['avg_distance_weight'] / count
#
#     print(f"  {edge_type}: {count}条边, 平均权重={avg_weight:.3f}, "
#           f"平均相关性={avg_corr:.3f}, 平均距离权重={avg_dist:.3f}")
#
# # 计算连通分量
# connected_components = list(nx.connected_components(G))
# print(f"连通分量数量: {len(connected_components)}")
# largest_component = max(connected_components, key=len)
# print(f"最大连通分量大小: {len(largest_component)}个节点")
#
# # 度分布
# degrees = [d for _, d in G.degree()]
# print(f"平均度: {np.mean(degrees):.2f}")
# print(f"最大度: {max(degrees)}")
#
# # ==================== 14. 保存结果（简化版） ====================
# print("\n=== 保存结果 ===")
#
# # 1. 保存PyG图（主要文件）
# torch.save(graph_data, '/datapool/home/info_wang/wanghui/file/graph/methylation_graph(5).pt')
# print("✓ PyG图已保存")
#
# # 2. 修复数据类型并保存NetworkX图（用于可视化分析）
# try:
#     import lxml
#
#     print("✓ lxml模块已安装")
# except ImportError:
#     import subprocess, sys
#
#     print("安装lxml模块...")
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "lxml", "--quiet"])
#     import lxml
#
#     print("✓ lxml模块安装完成")
#
# # 创建简化版的NetworkX图（只保存核心信息）
# G_simple = nx.Graph()
#
# # 添加节点
# for i, cpg in enumerate(cpg_list):
#     info = cpg_info_dict[cpg]
#     G_simple.add_node(i, cpg_id=str(cpg))  # 确保是字符串
#
# # 添加边（只保存权重）
# for i in range(edges[0].shape[0]):
#     source = edges[0][i]
#     target = edges[1][i]
#     if source != target:
#         weight = float(sparse_adj_matrix[source, target])
#         G_simple.add_edge(source, target, weight=weight)
#
# # 保存NetworkX图
# nx.write_graphml(G_simple, '/datapool/home/info_wang/wanghui/file/graph/methylation_graph(5).graphml')
# print("✓ NetworkX图已保存")
#
# # 3. 只保存一个汇总的CSV文件（包含所有重要信息）
# print("生成汇总文件...")
# summary_data = []
#
# # 节点信息
# for i, cpg in enumerate(cpg_list):
#     info = cpg_info_dict[cpg]
#     # 计算该节点的连接统计
#     neighbors = list(G_simple.neighbors(i))
#     if neighbors:
#         edge_weights = [G_simple[i][j]['weight'] for j in neighbors]
#         avg_edge_weight = float(np.mean(edge_weights))
#     else:
#         avg_edge_weight = 0.0
#
#     summary_data.append({
#         'node_id': i,
#         'cpg_id': cpg,
#         'site_type': str(info.site_type),
#         'chromosome': str(info.chromosome) if info.has_position else 'NA',
#         'position': int(info.position) if info.has_position else -1,
#         'regulatory_level': float(info.regulatory_weight),
#         'degree': len(neighbors),
#         'avg_edge_weight': avg_edge_weight,
#         'is_in_largest_component': i in largest_component
#     })
#
# # 创建DataFrame并保存
# summary_df = pd.DataFrame(summary_data)
# summary_df.to_csv('/datapool/home/info_wang/wanghui/file/graph/graph_summary(5).csv', index=False)
# print("✓ 汇总文件已保存")
#
# # 4. 输出关键统计信息到控制台
# print("\n" + "=" * 60)
# print("最终结果统计")
# print("=" * 60)
#
# print(f"\n📊 图结构信息:")
# print(f"   节点数: {graph_data.num_nodes}")
# print(f"   边数: {graph_data.num_edges}")
# print(f"   图密度: {nx.density(G_simple):.6f}")
# print(f"   平均度: {np.mean(degrees):.2f}")
# print(f"   最大连通分量: {len(largest_component)}个节点 ({len(largest_component) / n_cpgs * 100:.1f}%)")
#
# print(f"\n🔗 连接类型分布:")
# for edge_type, stats in type_connections.items():
#     count = stats['count']
#     percentage = count / G_simple.number_of_edges() * 100
#     print(f"   {edge_type}: {count}条边 ({percentage:.1f}%)")
#
# print(f"\n⚖️ 权重分布:")
# weight_values = combined_matrix[combined_matrix > 0]
# if len(weight_values) > 0:
#     weight_bins = np.histogram(weight_values, bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0])[0]
#     labels = ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']
#     total_edges = len(weight_values)
#     for label, count in zip(labels, weight_bins):
#         if count > 0:
#             print(f"   {label}: {count}条边 ({count / total_edges * 100:.1f}%)")
# else:
#     print("   没有有效的边权重")
#
# print(f"\n📈 各因子平均权重:")
# # 相关性权重（排除对角线）
# corr_values = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
# print(f"   平均相关性: {np.mean(corr_values):.4f}")
#
# # 距离权重（排除值为1的部分，即缺失位置的情况）
# dist_valid = dist_weight_matrix[dist_weight_matrix < 1.0]
# if len(dist_valid) > 0:
#     print(f"   平均距离权重: {np.mean(dist_valid):.4f}")
# else:
#     print(f"   平均距离权重: 无有效数据")
#
# # 调控权重
# reg_weights = [info.regulatory_weight for info in cpg_info_dict.values()]
# print(f"   平均调控权重: {np.mean(reg_weights):.4f}")
#
# # 综合权重
# if len(weight_values) > 0:
#     print(f"   平均综合权重: {np.mean(weight_values):.4f}")
#     print(f"   综合权重>0.5的边: {np.sum(weight_values > 0.5)}条")
#     print(f"   综合权重>0.8的边: {np.sum(weight_values > 0.8)}条")
#
# print(f"\n📍 文件保存位置:")
# print(f"   PyG图: /datapool/home/info_wang/wanghui/file/graph/methylation_graph(5).pt")
# print(f"   NetworkX图: /datapool/home/info_wang/wanghui/file/graph/methylation_graph(5).graphml")
# print(f"   汇总文件: /datapool/home/info_wang/wanghui/file/graph/graph_summary(5).csv")
#
# print("\n" + "=" * 60)
# print("✅ 分析完成！")
# print("=" * 60)
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import torch
from torch_geometric.data import Data
import networkx as nx
from tqdm import tqdm
import warnings
import os

warnings.filterwarnings('ignore')

# ==================== 1. 加载数据 ====================
print("正在加载数据...")
annotation_df = pd.read_csv("/datapool/home/info_wang/wanghui/file/450k_annotation.csv")
refseq_df = pd.read_csv("/datapool/home/info_wang/wanghui/file/ncbiRefSeq_hg19.csv")
df = pd.read_csv("/datapool/home/info_wang/wanghui/file/CPG/methylation_matrix_198.csv", index_col=0)

# 清理列名
df.columns = df.columns.str.strip()
methylation_data = df.drop(columns=['Age', 'Split', 'ID'], errors='ignore') # Added errors='ignore'
n_cpgs = len(methylation_data.columns)
n_samples = len(methylation_data)

print(f"甲基化数据: {n_samples}个样本, {n_cpgs}个CpG位点")

# ==================== 2. 整合注释信息 ====================
print("\n=== 整合注释信息 ===")
cpg_id_column = 'Name'

# 清理位点名称
annotation_df[cpg_id_column] = annotation_df[cpg_id_column].astype(str).str.strip()
available_cpgs = [cpg.strip() for cpg in methylation_data.columns]

# 筛选匹配的位点
cpg_positions = annotation_df[annotation_df[cpg_id_column].isin(available_cpgs)].copy()
print(f"在注释文件中找到 {len(cpg_positions)} 个匹配的CpG位点")

# 选择需要的列并重命名
cpg_positions = cpg_positions[['Name', 'chr', 'pos', 'Relation_to_Island', 'UCSC_RefGene_Name']].copy()
cpg_positions = cpg_positions.rename(columns={
    'Name': 'CpG_ID',
    'chr': 'CpG_Chromosome',
    'pos': 'Position',
    'Relation_to_Island': 'Region',
    'UCSC_RefGene_Name': 'Gene'
})

# ==================== 3. 处理基因区域信息 ====================
refseq_df.columns = refseq_df.columns.str.strip()
gene_regions = refseq_df[['name2', 'chrom', 'txStart', 'txEnd']].copy()
gene_regions = gene_regions.rename(columns={
    'name2': 'Gene',
    'chrom': 'Gene_Chromosome',
    'txStart': 'Start',
    'txEnd': 'End'
})

# 合并注释信息和基因区域信息
merged_df = pd.merge(cpg_positions, gene_regions, on='Gene', how='left')


# ==================== 4. 创建CpG位点信息类 ====================
class CpGSiteInfo:
    """存储每个CpG位点的完整信息"""

    def __init__(self, cpg_id):
        self.cpg_id = cpg_id
        self.has_position = False
        self.has_regulatory = False
        self.chromosome = None
        self.position = None
        self.regulatory_weight = None
        self.site_type = None  # 'A', 'B', 'C', 'D'

    def set_position_info(self, chromosome, position):
        """设置位置信息"""
        if pd.notna(chromosome) and pd.notna(position):
            self.chromosome = str(chromosome).strip()
            if not self.chromosome.startswith('chr'):
                self.chromosome = f'chr{self.chromosome}'
            self.position = int(position)
            self.has_position = True

    def set_regulatory_info(self, region, gene_info):
        """根据注释信息设置调控权重"""
        # 默认权重
        self.regulatory_weight = 0.55  # 保守中间值
        self.has_regulatory = False

        # 如果有区域信息
        if pd.notna(region):
            region_str = str(region).upper()

            # 第1级：核心调控区域 (权重=1.0)
            if any(keyword in region_str for keyword in ['ISLAND', 'TSST1500']):
                self.regulatory_weight = 1.0
                self.has_regulatory = True

            # 第2级：潜在调控区域 (权重=0.7)
            elif any(keyword in region_str for keyword in ['SHORE', 'SHELF']):
                self.regulatory_weight = 0.7
                self.has_regulatory = True

            # 第3级：非调控区域 (权重=0.4)
            elif 'OPENSEA' in region_str:
                self.regulatory_weight = 0.4
                self.has_regulatory = True

        # 如果有基因信息，可以提升权重
        if pd.notna(gene_info) and self.has_regulatory:
            # 如果有基因关联，稍微提升权重
            self.regulatory_weight = min(1.0, self.regulatory_weight + 0.1)

    def determine_site_type(self):
        """确定位点类型"""
        if self.has_position and self.has_regulatory:
            self.site_type = 'A'
        elif not self.has_position and self.has_regulatory:
            self.site_type = 'B'
        elif self.has_position and not self.has_regulatory:
            self.site_type = 'C'
        else:  # not has_position and not has_regulatory
            self.site_type = 'D'
        return self.site_type


# ==================== 5. 初始化所有CpG位点信息 ====================
print("\n=== 初始化CpG位点信息 ===")
cpg_info_dict = {}

for cpg in methylation_data.columns:
    info = CpGSiteInfo(cpg)

    # 查找位置信息
    cpg_data = merged_df[merged_df['CpG_ID'] == cpg]
    if not cpg_data.empty:
        info.set_position_info(
            cpg_data.iloc[0]['CpG_Chromosome'],
            cpg_data.iloc[0]['Position']
        )

        # 设置调控信息
        info.set_regulatory_info(
            cpg_data.iloc[0]['Region'],
            cpg_data.iloc[0]['Gene']
        )

    # 确定位点类型
    info.determine_site_type()
    cpg_info_dict[cpg] = info

# 统计位点类型
type_counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
for info in cpg_info_dict.values():
    type_counts[info.site_type] += 1

print(f"位点类型统计:")
print(f"  A类(有位置+有调控): {type_counts['A']} ({type_counts['A'] / n_cpgs * 100:.1f}%)")
print(f"  B类(无位置+有调控): {type_counts['B']} ({type_counts['B'] / n_cpgs * 100:.1f}%)")
print(f"  C类(有位置+无调控): {type_counts['C']} ({type_counts['C'] / n_cpgs * 100:.1f}%)")
print(f"  D类(无位置+无调控): {type_counts['D']} ({type_counts['D'] / n_cpgs * 100:.1f}%)")

# ==================== 6. 计算相关性权重矩阵 ====================
print("\n=== 计算相关性权重矩阵 ===")


def calculate_correlation_matrix_optimized(data_matrix):
    """优化版本的相关性计算"""
    n_features = data_matrix.shape[1]
    corr_matrix = np.ones((n_features, n_features))

    # 转换为numpy数组并转置为(特征数, 样本数)
    data_array = data_matrix.values.T

    print("计算皮尔森相关系数...")
    for i in tqdm(range(n_features), desc="处理位点"):
        for j in range(i + 1, n_features):
            # 获取两个位点的数据
            x = data_array[i]
            y = data_array[j]

            # 去除缺失值
            mask = ~(np.isnan(x) | np.isnan(y))
            x_valid = x[mask]
            y_valid = y[mask]

            # 至少需要3个有效样本
            if len(x_valid) < 3:
                corr = 0.0
            else:
                # 计算皮尔森相关系数
                corr, _ = pearsonr(x_valid, y_valid)
                if np.isnan(corr):
                    corr = 0.0

            # 存储到矩阵（使用绝对值）
            corr_matrix[i, j] = abs(corr)
            corr_matrix[j, i] = abs(corr)

    return corr_matrix


# 计算相关性矩阵
corr_matrix = calculate_correlation_matrix_optimized(methylation_data)
print(f"相关性矩阵计算完成，形状: {corr_matrix.shape}")

# ==================== 7. 计算距离权重矩阵 ====================
print("\n=== 计算距离权重矩阵 ===")


def calculate_distance_weight(distance):
    """分段计算距离权重"""
    if distance < 1000:  # < 1kb
        return np.exp(-distance / 500)
    elif distance < 10000:  # 1kb - 10kb
        return 0.9 * np.exp(-(distance - 1000) / 5000)
    elif distance < 100000:  # 10kb - 100kb
        return 0.6 * np.exp(-(distance - 10000) / 30000)
    elif distance < 1000000:  # 100kb - 1Mb
        return 0.3 * np.exp(-(distance - 100000) / 300000)
    else:  # > 1Mb
        return 0.1 * np.exp(-(distance - 1000000) / 1000000)


def calculate_distance_weight_matrix(cpg_info_dict, cpg_list):
    """计算距离权重矩阵"""
    n = len(cpg_list)
    dist_weight_matrix = np.ones((n, n))  # 默认值为1（中性）

    # 创建位点索引映射
    cpg_to_idx = {cpg: i for i, cpg in enumerate(cpg_list)}

    # 只对有位置信息的位点进行计算
    valid_cpgs = [cpg for cpg in cpg_list if cpg_info_dict[cpg].has_position]
    print(f"有位置信息的位点数量: {len(valid_cpgs)}/{n}")

    # 计算同染色体位点的距离
    chromosome_groups = {}
    for cpg in valid_cpgs:
        info = cpg_info_dict[cpg]
        chrom = info.chromosome
        if chrom not in chromosome_groups:
            chromosome_groups[chrom] = []
        chromosome_groups[chrom].append((cpg, info.position, cpg_to_idx[cpg]))

    print(f"分组到 {len(chromosome_groups)} 条染色体")

    # 对每条染色体内部的位点计算距离
    for chrom, sites in tqdm(chromosome_groups.items(), desc="计算染色体内部距离"):
        if len(sites) < 2:
            continue

        # 按位置排序
        sites.sort(key=lambda x: x[1])

        for i in range(len(sites)):
            cpg1, pos1, idx1 = sites[i]
            for j in range(i + 1, len(sites)):
                cpg2, pos2, idx2 = sites[j]
                distance = abs(pos1 - pos2)
                weight = calculate_distance_weight(distance)
                dist_weight_matrix[idx1, idx2] = weight
                dist_weight_matrix[idx2, idx1] = weight

    # 不同染色体的位点对，权重设为0.001
    different_chromosome_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            info_i = cpg_info_dict[cpg_list[i]]
            info_j = cpg_info_dict[cpg_list[j]]

            if info_i.has_position and info_j.has_position:
                if info_i.chromosome != info_j.chromosome:
                    dist_weight_matrix[i, j] = 0.001
                    dist_weight_matrix[j, i] = 0.001
                    different_chromosome_pairs += 1

    print(f"不同染色体位点对数量: {different_chromosome_pairs}")
    return dist_weight_matrix


# 计算距离权重矩阵
cpg_list = list(methylation_data.columns)
dist_weight_matrix = calculate_distance_weight_matrix(cpg_info_dict, cpg_list)

# ==================== 8. 计算调控权重矩阵 ====================
print("\n=== 计算调控权重矩阵 ===")


def calculate_regulatory_weight_matrix(cpg_info_dict, cpg_list):
    """计算调控权重矩阵"""
    n = len(cpg_list)
    reg_weight_matrix = np.ones((n, n))

    # 获取每个位点的调控权重
    reg_weights = np.array([cpg_info_dict[cpg].regulatory_weight for cpg in cpg_list])

    # 计算调控权重矩阵：使用几何平均 R_ij = sqrt(r_i * r_j)
    for i in range(n):
        for j in range(n):
            if i == j:
                reg_weight_matrix[i, j] = reg_weights[i]
            else:
                reg_weight_matrix[i, j] = np.sqrt(reg_weights[i] * reg_weights[j])

    return reg_weight_matrix


# 计算调控权重矩阵
reg_weight_matrix = calculate_regulatory_weight_matrix(cpg_info_dict, cpg_list)

# ==================== 9. 计算综合权重矩阵 ====================
print("\n=== 计算综合权重矩阵 ===")


def calculate_combined_weight(corr_weight, dist_weight, reg_weight,
                              info1, info2):
    """根据信息完整性动态计算综合权重"""

    # 确定有效因子数
    n_factors = 1  # 相关性总是有
    factors = [corr_weight]

    # 距离因子
    if info1.has_position and info2.has_position:
        n_factors += 1
        factors.append(dist_weight)
    else:
        # 无位置信息时，距离权重设为1（中性）
        factors.append(1.0)

    # 调控因子
    if info1.has_regulatory or info2.has_regulatory:
        n_factors += 1
        factors.append(reg_weight)
    else:
        # 两者都无调控信息时，调控权重设为0.55（保守值）
        factors.append(0.55)

    # 计算几何平均
    if n_factors == 3:
        return np.prod(factors) ** (1 / 3)
    elif n_factors == 2:
        # 只使用相关性和调控（或无位置但有调控的情况）
        return (corr_weight * factors[2]) ** (1 / 2)
    else:  # n_factors == 1
        return corr_weight


def calculate_combined_weight_matrix(corr_matrix, dist_matrix, reg_matrix,
                                     cpg_info_dict, cpg_list):
    """计算综合权重矩阵"""
    n = len(cpg_list)
    combined_matrix = np.zeros((n, n))
    info_completeness = np.zeros((n, n), dtype='U3')  # 存储信息完整性

    print("计算综合权重...")
    for i in tqdm(range(n), desc="处理位点"):
        info_i = cpg_info_dict[cpg_list[i]]
        for j in range(i, n):
            info_j = cpg_info_dict[cpg_list[j]]

            # 计算综合权重
            weight = calculate_combined_weight(
                corr_matrix[i, j],
                dist_matrix[i, j],
                reg_matrix[i, j],
                info_i,
                info_j
            )

            combined_matrix[i, j] = weight
            combined_matrix[j, i] = weight

            # 记录信息完整性
            completeness = f"{info_i.site_type}{info_j.site_type}"
            info_completeness[i, j] = completeness
            info_completeness[j, i] = completeness

    return combined_matrix, info_completeness


# 计算综合权重矩阵
combined_matrix, info_completeness = calculate_combined_weight_matrix(
    corr_matrix, dist_weight_matrix, reg_weight_matrix,
    cpg_info_dict, cpg_list
)

# ==================== 10. 稀疏化处理 ====================
print("\n=== 稀疏化处理 ===")


def sparse_top_k(matrix, k=20):
    """每个节点只保留权重最高的k个连接"""
    n = matrix.shape[0]
    sparse_matrix = np.zeros_like(matrix)

    for i in range(n):
        # 获取当前节点的所有权重
        weights = matrix[i].copy()
        weights[i] = -np.inf  # 排除自身

        # 找出top-k个邻居
        top_k_indices = np.argsort(weights)[-k:]
        sparse_matrix[i, top_k_indices] = matrix[i, top_k_indices]

    # 确保对称性
    sparse_matrix = np.maximum(sparse_matrix, sparse_matrix.T)

    # 计算稀疏化后的统计
    edge_count = np.sum(sparse_matrix > 0) // 2
    print(f"稀疏化后边数量: {edge_count}")
    print(f"平均每个节点连接数: {2 * edge_count / n:.2f}")

    return sparse_matrix


# 稀疏化处理
k_neighbors = 20
sparse_adj_matrix = sparse_top_k(combined_matrix, k=k_neighbors)

# ==================== 11. 构建图结构 ====================
print("\n=== 构建图结构 ===")

# 提取边
edges = np.nonzero(sparse_adj_matrix)
edge_index = torch.tensor(edges, dtype=torch.long)

# 节点特征：甲基化数据转置
node_features = torch.tensor(methylation_data.values.T, dtype=torch.float)

# 边权重
edge_weights = torch.tensor(sparse_adj_matrix[edges], dtype=torch.float)

# 创建PyG图数据
graph_data = Data(
    x=node_features,
    edge_index=edge_index,
    edge_attr=edge_weights.unsqueeze(1)  # 转换为列向量
)

print(f"图数据信息:")
print(f"- 节点数: {graph_data.num_nodes}")
print(f"- 边数: {graph_data.num_edges}")
print(f"- 节点特征维度: {graph_data.num_node_features}")
print(f"- 边特征维度: {graph_data.num_edge_features}")

# ==================== 12. 创建NetworkX图（用于分析） ====================
print("\n=== 创建NetworkX图 ===")
G = nx.Graph()

# 添加节点
for i, cpg in enumerate(cpg_list):
    info = cpg_info_dict[cpg]
    G.add_node(i,
               cpg_id=cpg,
               site_type=info.site_type,
               chromosome=info.chromosome if info.has_position else None,
               position=info.position if info.has_position else None,
               regulatory_weight=info.regulatory_weight)

# 添加边
for i in range(edges[0].shape[0]):
    source = edges[0][i]
    target = edges[1][i]
    if source != target:
        weight = sparse_adj_matrix[source, target]
        # 获取信息完整性
        completeness = info_completeness[source, target]
        G.add_edge(source, target,
                   weight=weight,
                   info_completeness=completeness,
                   correlation=corr_matrix[source, target],
                   distance_weight=dist_weight_matrix[source, target],
                   regulatory_weight=reg_weight_matrix[source, target])

print(f"NetworkX图信息:")
print(f"- 节点数: {G.number_of_nodes()}")
print(f"- 边数: {G.number_of_edges()}")
print(f"- 图密度: {nx.density(G):.6f}")

# ==================== 13. 图统计分析 ====================
print("\n=== 图统计分析 ===")

# 按位点类型统计连接
type_connections = {}
for u, v, data in G.edges(data=True):
    u_type = G.nodes[u]['site_type']
    v_type = G.nodes[v]['site_type']
    edge_type = f"{u_type}-{v_type}"

    if edge_type not in type_connections:
        type_connections[edge_type] = {
            'count': 0,
            'total_weight': 0,
            'avg_correlation': 0,
            'avg_distance_weight': 0
        }

    type_connections[edge_type]['count'] += 1
    type_connections[edge_type]['total_weight'] += data['weight']
    type_connections[edge_type]['avg_correlation'] += data['correlation']
    type_connections[edge_type]['avg_distance_weight'] += data['distance_weight']

print("按连接类型统计:")
for edge_type, stats in sorted(type_connections.items()):
    count = stats['count']
    avg_weight = stats['total_weight'] / count
    avg_corr = stats['avg_correlation'] / count
    avg_dist = stats['avg_distance_weight'] / count

    print(f"  {edge_type}: {count}条边, 平均权重={avg_weight:.3f}, "
          f"平均相关性={avg_corr:.3f}, 平均距离权重={avg_dist:.3f}")

# 计算连通分量
connected_components = list(nx.connected_components(G))
print(f"连通分量数量: {len(connected_components)}")
largest_component = max(connected_components, key=len)
print(f"最大连通分量大小: {len(largest_component)}个节点")

# 度分布
degrees = [d for _, d in G.degree()]
print(f"平均度: {np.mean(degrees):.2f}")
print(f"最大度: {max(degrees)}")
# ==================== 14. 保存结果（含全局图指标） ====================
print("\n=== 保存结果 ===")

# 1. 保存PyG图（主要文件）
torch.save(graph_data, '/datapool/home/info_wang/wanghui/file/graph/methylation_graph(20).pt')
print("✓ PyG图已保存")

# 2. 保存NetworkX图（用于可视化分析）
# 创建简化版NetworkX图
G_simple = nx.Graph()
for i, cpg in enumerate(cpg_list):
    info = cpg_info_dict[cpg]
    G_simple.add_node(i, cpg_id=str(cpg))
for i in range(edges[0].shape[0]):
    source = edges[0][i]
    target = edges[1][i]
    if source != target:
        weight = float(sparse_adj_matrix[source, target])
        G_simple.add_edge(source, target, weight=weight)

nx.write_graphml(G_simple, '/datapool/home/info_wang/wanghui/file/graph/methylation_graph(20).graphml')
print("✓ NetworkX图已保存")

# 3. 生成增强版汇总文件（包含节点信息 + 全局图指标）
print("正在计算图指标并生成汇总文件...")

# --- A. 计算节点级指标 ---
print("  - 计算节点聚类系数...")
clustering_coeffs = nx.clustering(G_simple, weight='weight')

print("  - 计算节点中心性...")
try:
    eigenvector_centrality = nx.eigenvector_centrality(G_simple, weight='weight', max_iter=1000)
except:
    eigenvector_centrality = nx.degree_centrality(G_simple)

# --- B. 计算全局图指标 (新增功能) ---
print("  - 计算全局图指标...")
global_num_nodes = G_simple.number_of_nodes()
global_num_edges = G_simple.number_of_edges()
global_density = nx.density(G_simple)
global_avg_clustering = nx.average_clustering(G_simple, weight='weight')
global_avg_degree = 2 * global_num_edges / global_num_nodes if global_num_nodes > 0 else 0

# 计算最大连通分量占比
connected_components = list(nx.connected_components(G_simple))
largest_component = max(connected_components, key=len)
global_largest_component_ratio = len(largest_component) / global_num_nodes

summary_data = []

# --- C. 整合数据 ---
for i, cpg in enumerate(cpg_list):
    info = cpg_info_dict[cpg]
    neighbors = list(G_simple.neighbors(i))

    if neighbors:
        edge_weights = [G_simple[i][j]['weight'] for j in neighbors]
        avg_edge_weight = float(np.mean(edge_weights))
    else:
        avg_edge_weight = 0.0

    # 构造每一行数据
    row_data = {
        # --- 节点自身信息 ---
        'node_id': i,
        'cpg_id': cpg,
        'site_type': str(info.site_type),
        'chromosome': str(info.chromosome) if info.has_position else 'NA',
        'position': int(info.position) if info.has_position else -1,
        'regulatory_level': float(info.regulatory_weight),

        # --- 节点拓扑特征 ---
        'degree': len(neighbors),
        'avg_edge_weight': avg_edge_weight,
        'clustering_coefficient': clustering_coeffs.get(i, 0.0),
        'eigenvector_centrality': eigenvector_centrality.get(i, 0.0),
        'is_in_largest_component': i in largest_component,

        # --- 全局图特征 (每行重复，方便读取) ---
        'graph_total_nodes': global_num_nodes,
        'graph_total_edges': global_num_edges,
        'graph_density': global_density,
        'graph_avg_degree': global_avg_degree,
        'graph_avg_clustering': global_avg_clustering,
        'graph_largest_comp_ratio': global_largest_component_ratio
    }
    summary_data.append(row_data)

# 创建并保存
summary_df = pd.DataFrame(summary_data)
output_csv = '/datapool/home/info_wang/wanghui/file/graph/graph_summary(20).csv'
summary_df.to_csv(output_csv, index=False)
print(f"✓ 汇总文件已保存: {output_csv}")

# 4. 控制台打印确认
print("\n" + "=" * 60)
print("📊 新增的全局指标预览 (已写入CSV):")
print(f"   - 总节点数: {global_num_nodes}")
print(f"   - 总边数:   {global_num_edges}")
print(f"   - 图密度:   {global_density:.4f}")
print(f"   - 平均度:   {global_avg_degree:.2f}")
print(f"   - 平均聚类: {global_avg_clustering:.4f}")
print("=" * 60)