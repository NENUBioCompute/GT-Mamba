import pandas as pd
import numpy as np
import os
import json
from pathlib import Path
from sklearn.linear_model import LinearRegression
import warnings

warnings.filterwarnings('ignore')


def check_health_status(pheno_df, dataset_name):
    """
    检查样本健康状态
    """
    health_indicators = {
        'healthy_samples': [],
        'unhealthy_samples': [],
        'unknown_samples': [],
        'health_criteria_used': []
    }

    all_samples = pheno_df.index.tolist()

    # 健康状态列名，包含 Disease
    health_columns = [
        'Disease', 'disease', 'health_status', 'health', 'disease_status',
        'condition', 'diagnosis', 'phenotype', 'status', 'group',
        'sample_type', 'type', 'category', 'clinical_status'
    ]

    available_health_columns = [col for col in health_columns if col in pheno_df.columns]

    if not available_health_columns:
        print(f"  警告: 数据集 {dataset_name} 没有找到健康状态列，默认所有样本为健康")
        health_indicators['healthy_samples'] = all_samples
        health_indicators['health_criteria_used'] = ['no_health_info_default_healthy']
        return health_indicators

    health_indicators['health_criteria_used'] = available_health_columns

    print(f"  找到健康状态列: {available_health_columns}")

    # 显示该列的唯一值，帮助调试
    main_health_col = available_health_columns[0]
    unique_values = pheno_df[main_health_col].unique()
    print(f"  '{main_health_col}' 列的唯一值: {list(unique_values)}")

    # 健康关键词（不区分大小写）
    healthy_keywords = [
        'healthy', 'normal', 'control', 'ctrl', 'no disease',
        'nondisease', 'none', 'naive', 'wild type', 'wt', 'healthy control',
        'normal control', 'reference', 'baseline', '0', 'no', 'false',
        'negative', 'non-disease', 'health'
    ]

    # 疾病关键词
    disease_keywords = [
        'disease', 'cancer', 'tumor', 'alzheimer', 'diabetes',
        'parkinson', 'sick', 'case', 'patient', 'ad', 'cvd',
        'cardiovascular', 'carcinoma', 'leukemia', 'lymphoma',
        'melanoma', 'adenocarcinoma', 'tumor', 'malignant',
        '1', 'yes', 'true', 'positive', 'affected'
    ]

    for sample in all_samples:
        sample_health_info = []

        for health_col in available_health_columns:
            health_value = str(pheno_df.loc[sample, health_col]).lower().strip()
            sample_health_info.append(health_value)

        # 判断健康状态
        health_info_str = ' '.join(sample_health_info)

        is_healthy = any(keyword in health_info_str for keyword in healthy_keywords)
        is_unhealthy = any(keyword in health_info_str for keyword in disease_keywords)

        if is_healthy and not is_unhealthy:
            health_indicators['healthy_samples'].append(sample)
        elif is_unhealthy:
            health_indicators['unhealthy_samples'].append(sample)
        else:
            health_indicators['unknown_samples'].append(sample)

    print(f"  健康样本: {len(health_indicators['healthy_samples'])}")
    print(f"  疾病样本: {len(health_indicators['unhealthy_samples'])}")
    print(f"  未知状态: {len(health_indicators['unknown_samples'])}")

    return health_indicators


def load_flowsorted_reference():
    """
    加载 FlowSorted.Blood.450k 参考数据集
    """
    reference_file = '/mnt/sfs_turbo/wanghui/protein/file/FlowSorted.Blood.450k.csv'
    print(f"从 {reference_file} 加载参考数据集...")

    # 检查文件是否存在
    if not os.path.exists(reference_file):
        print(f"❌ 参考数据文件不存在: {reference_file}")
        return None, None

    try:
        # 加载参考数据
        print(f"正在加载参考数据文件...")
        reference_data = pd.read_csv(reference_file, index_col=0)
        print(f"✓ 参考数据加载成功: {reference_data.shape}")

        # 检查数据格式
        print(f"数据前5行5列:")
        print(reference_data.iloc[:5, :5])

        # 判断这是否是参考矩阵格式
        # FlowSorted.Blood.450k 通常是样本在列，探针在行
        # 列名应该包含细胞类型信息

        # 创建表型数据（从列名推断细胞类型）
        reference_pheno = pd.DataFrame(index=reference_data.columns)

        # 从列名推断细胞类型
        cell_types = []
        for col in reference_data.columns:
            col_lower = col.lower()
            if 'cd8' in col_lower or 'tcytotoxic' in col_lower:
                cell_types.append('CD8T')
            elif 'cd4' in col_lower or 'thelper' in col_lower:
                cell_types.append('CD4T')
            elif 'nk' in col_lower:
                cell_types.append('NK')
            elif 'bcell' in col_lower or 'b_' in col_lower or 'b.lymphocyte' in col_lower:
                cell_types.append('Bcell')
            elif 'mono' in col_lower:
                cell_types.append('Mono')
            elif 'gran' in col_lower or 'neut' in col_lower:
                cell_types.append('Granulocytes')
            elif 'eos' in col_lower:
                cell_types.append('Eosinophils')
            elif 'baso' in col_lower:
                cell_types.append('Basophils')
            else:
                # 尝试其他常见细胞类型
                for ct in ['Neutrophil', 'Lymphocyte', 'Monocyte', 'Eosinophil', 'Basophil']:
                    if ct.lower() in col_lower:
                        cell_types.append(ct)
                        break
                else:
                    cell_types.append('Unknown')

        reference_pheno['CellType'] = cell_types
        print(f"推断的细胞类型: {list(set(cell_types))}")

        return reference_data, reference_pheno

    except Exception as e:
        print(f"❌ 加载参考数据失败: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def load_beta_data(data_dir, dataset_name):
    """
    加载beta数据，支持分片文件
    """
    print(f"  加载 {dataset_name} 的beta数据...")

    # 检查是否存在分片文件
    beta_files = []

    # 首先检查是否存在完整的beta文件
    full_beta_file = os.path.join(data_dir, f"{dataset_name}_beta.csv")
    if os.path.exists(full_beta_file):
        print(f"    找到完整文件: {full_beta_file}")
        try:
            beta_df = pd.read_csv(full_beta_file, index_col=0)
            print(f"    ✓ Beta数据加载成功: {beta_df.shape}")
            return beta_df
        except Exception as e:
            print(f"    ❌ 加载完整文件失败: {e}")

    # 检查分片文件
    for i in range(1, 20):  # 假设最多20个分片
        beta_file = os.path.join(data_dir, f"{dataset_name}_beta_{i}.csv")
        if os.path.exists(beta_file):
            beta_files.append(beta_file)
            print(f"    找到分片文件: {beta_file}")

    if not beta_files:
        print(f"    ❌ 未找到 {dataset_name} 的beta文件")
        return None

    # 加载并合并所有分片
    beta_dfs = []
    for beta_file in beta_files:
        try:
            beta_part = pd.read_csv(beta_file, index_col=0)
            beta_dfs.append(beta_part)
            print(f"    加载分片: {beta_part.shape}")
        except Exception as e:
            print(f"    ❌ 加载分片 {beta_file} 失败: {e}")

    if not beta_dfs:
        print(f"    ❌ 所有分片加载失败")
        return None

    # 合并所有分片
    if len(beta_dfs) == 1:
        beta_df = beta_dfs[0]
    else:
        # 按列合并所有分片
        try:
            beta_df = pd.concat(beta_dfs, axis=1)
            print(f"    ✓ 合并成功，最终形状: {beta_df.shape}")
        except Exception as e:
            print(f"    ❌ 合并分片失败: {e}")
            return None

    print(f"    ✓ 最终beta数据形状: {beta_df.shape}")
    return beta_df


def calculate_cell_proportions_with_reference(beta_df, dataset_name, reference_beta, reference_pheno):
    """
    使用参考数据集计算细胞比例
    """
    print("  使用 FlowSorted.Blood.450k 计算细胞比例...")

    try:
        # 1. 找到共有的CpG位点
        common_probes = beta_df.index.intersection(reference_beta.index)
        print(f"  共有CpG位点: {len(common_probes)} 个")

        if len(common_probes) < 100:
            print("  ⚠️ 共有CpG位点太少，细胞比例计算可能不准确")
            return None

        # 2. 提取共有的数据
        beta_common = beta_df.loc[common_probes]
        reference_common = reference_beta.loc[common_probes]

        # 3. 获取参考数据集的细胞类型信息
        if 'CellType' in reference_pheno.columns:
            cell_types = reference_pheno['CellType'].unique()
        else:
            # 使用常见的血液细胞类型
            cell_types = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Granulocytes']

        print(f"  参考细胞类型: {list(cell_types)}")

        # 4. 为每种细胞类型创建平均参考谱
        reference_means = {}
        valid_cell_types = []

        for cell_type in cell_types:
            # 找到属于该细胞类型的样本
            cell_samples = []

            if 'CellType' in reference_pheno.columns:
                cell_samples = reference_pheno[reference_pheno['CellType'] == cell_type].index
            else:
                # 从列名推断
                cell_samples = [col for col in reference_common.columns if str(cell_type).lower() in col.lower()]

            if len(cell_samples) > 0:
                cell_data = reference_common[cell_samples]
                # 移除缺失值过多的位点
                valid_probes = cell_data.isnull().mean(axis=1) < 0.5
                if valid_probes.sum() > 100:
                    reference_means[cell_type] = cell_data.loc[valid_probes].mean(axis=1)
                    valid_cell_types.append(cell_type)

        if not valid_cell_types:
            print("  ⚠️ 没有找到有效的细胞类型参考谱，使用默认细胞类型")
            valid_cell_types = ['CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Granulocytes']
            # 创建简单的参考谱
            for cell_type in valid_cell_types:
                reference_means[cell_type] = reference_common.mean(axis=1)

        print(f"  有效细胞类型: {valid_cell_types}")

        # 5. 创建参考矩阵
        reference_matrix = pd.DataFrame(reference_means)
        print(f"  参考矩阵形状: {reference_matrix.shape}")

        # 6. 处理缺失值
        reference_matrix = reference_matrix.fillna(reference_matrix.mean())

        # 7. 使用线性回归进行去卷积
        cell_proportions = pd.DataFrame(index=beta_common.columns, columns=valid_cell_types)

        samples_processed = 0
        for i, sample in enumerate(beta_common.columns):
            if i % 50 == 0:  # 每50个样本打印一次进度
                print(f"    正在处理样本 {i + 1}/{len(beta_common.columns)}")

            # 对于每个样本，使用线性回归估计细胞比例
            X = reference_matrix.values
            y = beta_common[sample].values

            # 处理y中的缺失值
            valid_indices = ~np.isnan(y)
            if valid_indices.sum() < 100:
                continue

            X_valid = X[valid_indices]
            y_valid = y[valid_indices]

            # 使用非负线性回归
            try:
                model = LinearRegression(fit_intercept=False, positive=True)
                model.fit(X_valid, y_valid)

                proportions = model.coef_
                proportions = np.clip(proportions, 0, 1)
                if proportions.sum() > 0:
                    proportions = proportions / proportions.sum()

                cell_proportions.loc[sample] = proportions
                samples_processed += 1
            except:
                # 如果回归失败，使用均匀分布
                cell_proportions.loc[sample] = [1.0 / len(valid_cell_types)] * len(valid_cell_types)
                samples_processed += 1

        print(f"  ✓ 成功处理 {samples_processed}/{len(beta_common.columns)} 个样本")
        return cell_proportions.T

    except Exception as e:
        print(f"  ✗ 细胞比例计算失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def process_single_dataset_with_cell_proportions(dataset, data_dir, pheno_output_dir, reference_beta, reference_pheno):
    """
    处理单个数据集并计算细胞比例
    """
    print(f"\n{'=' * 50}")
    print(f"处理数据集: {dataset}")
    print(f"{'=' * 50}")

    try:
        # 读取表型数据
        pheno_file = os.path.join(data_dir, f"{dataset}_pheno.csv")
        if not os.path.exists(pheno_file):
            print(f"  ❌ 表型文件不存在: {pheno_file}")
            return {'status': 'error', 'error_message': f'表型文件不存在: {pheno_file}'}

        pheno_df = pd.read_csv(pheno_file)
        print(f"  原始表型数据: {pheno_df.shape}")

        # 检查表型数据的索引列
        if 'ID' in pheno_df.columns:
            # 使用ID列作为索引
            pheno_df = pheno_df.set_index('ID')
            print(f"  使用ID列作为索引，设置后形状: {pheno_df.shape}")
        elif 'Sample' in pheno_df.columns:
            # 使用Sample列作为索引
            pheno_df = pheno_df.set_index('Sample')
            print(f"  使用Sample列作为索引，设置后形状: {pheno_df.shape}")
        else:
            # 使用第一列作为索引
            first_col = pheno_df.columns[0]
            if 'Unnamed' in first_col:
                pheno_df = pheno_df.set_index(pheno_df.columns[0])
                print(f"  使用第一列作为索引，设置后形状: {pheno_df.shape}")
            else:
                print(f"  使用默认索引，形状: {pheno_df.shape}")

        # 加载beta数据
        beta_df = load_beta_data(data_dir, dataset)
        if beta_df is None:
            print(f"  ❌ 无法加载beta数据，跳过 {dataset}")
            return {'status': 'error', 'error_message': '无法加载beta数据'}

        print(f"  Beta数据: {beta_df.shape}")

        # 数据对齐 - 详细检查
        print("  数据对齐检查...")
        pheno_samples = set(pheno_df.index)
        beta_samples = set(beta_df.columns)

        print(f"  表型数据样本数: {len(pheno_samples)}")
        print(f"  Beta数据样本数: {len(beta_samples)}")

        common_samples = pheno_samples & beta_samples
        print(f"  共同样本数量: {len(common_samples)}")

        if len(common_samples) == 0:
            print(f"  ❌ 错误: 没有共同的样本")
            # 显示一些样本对比
            print(f"  表型数据前5个样本: {list(pheno_samples)[:5]}")
            print(f"  Beta数据前5个样本: {list(beta_samples)[:5]}")
            return {'status': 'error', 'error_message': '没有共同的样本'}

        beta_df_aligned = beta_df[list(common_samples)]
        pheno_df_aligned = pheno_df.loc[list(common_samples)]

        print(f"  对齐后 - 表型: {pheno_df_aligned.shape}, 甲基化: {beta_df_aligned.shape}")

        # 健康样本检测
        print("  步骤3: 健康样本检测...")
        health_indicators = check_health_status(pheno_df_aligned, dataset)

        # 只保留健康样本
        healthy_samples = health_indicators['healthy_samples']
        if not healthy_samples:
            print(f"  ⚠️ 警告: 数据集 {dataset} 没有健康样本，跳过")
            return {'status': 'error', 'error_message': '没有健康样本'}

        pheno_healthy = pheno_df_aligned.loc[healthy_samples]
        beta_healthy = beta_df_aligned[healthy_samples]

        print(f"  健康样本 - 表型: {pheno_healthy.shape}, 甲基化: {beta_healthy.shape}")

        # 计算细胞比例（如果参考数据可用）
        cell_proportions = None
        cell_types = []
        if reference_beta is not None:
            print("  步骤4: 计算细胞比例...")
            cell_proportions = calculate_cell_proportions_with_reference(
                beta_healthy, dataset, reference_beta, reference_pheno
            )

            if cell_proportions is not None:
                for cell_type in cell_proportions.index:
                    pheno_healthy[f'cell_{cell_type}'] = cell_proportions.loc[cell_type].values
                cell_types = cell_proportions.index.tolist()
                print(f"  ✓ 细胞比例添加成功")
            else:
                print(f"  ⚠️ 细胞比例计算失败")
        else:
            print("  ⚠️ 无参考数据，跳过细胞比例计算")

        # 添加元数据
        pheno_healthy['health_status'] = 'healthy'
        pheno_healthy['dataset'] = dataset
        pheno_healthy['health_check_method'] = str(health_indicators['health_criteria_used'])
        pheno_healthy['cell_proportion_calculated'] = cell_proportions is not None
        pheno_healthy['processing_timestamp'] = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')

        # 保存新的表型文件
        pheno_output_path = os.path.join(pheno_output_dir, f"{dataset}_pheno_healthy_with_celltypes.csv")
        pheno_healthy.to_csv(pheno_output_path)
        print(f"  ✓ 保存表型文件: {pheno_output_path}")

        # 记录处理结果
        return {
            'status': 'success',
            'pheno_file': pheno_output_path,
            'samples': healthy_samples,
            'sample_count': len(healthy_samples),
            'cell_types': cell_types,
            'has_cell_proportions': cell_proportions is not None,
            'health_stats': {
                'total_samples': len(common_samples),
                'healthy_samples': len(healthy_samples),
                'unhealthy_samples': len(health_indicators['unhealthy_samples']),
                'unknown_samples': len(health_indicators['unknown_samples']),
                'health_columns_used': health_indicators['health_criteria_used']
            }
        }

    except Exception as e:
        print(f"  ❌ 处理数据集 {dataset} 时出错: {e}")
        import traceback
        traceback.print_exc()
        return {
            'status': 'error',
            'error_message': str(e)
        }


def process_all_datasets_with_cell_proportions():
    """
    处理所有数据集并计算细胞比例，支持分片文件
    """
    # 配置路径
    data_dir = '/mnt/sfs_turbo/wanghui/protein/data/GSE-data'
    pheno_output_dir = '/mnt/sfs_turbo/wanghui/protein/data/GSE-data'
    other_files_dir = '/mnt/sfs_turbo/wanghui/protein/file'

    # 创建输出目录
    Path(pheno_output_dir).mkdir(exist_ok=True)
    Path(other_files_dir).mkdir(exist_ok=True)

    # 所有目标数据集
    target_datasets = [
        'GSE40279', 'GSE72777', 'GSE77445', 'GSE61496',
        'GSE50660', 'GSE53128', 'GSE55763',
        'GSE60132', 'GSE64495', 'GSE65638',  'GSE72773',
        'GSE72775', 'GSE73103', 'GSE87571'
    ]

    print("=" * 70)
    print("开始处理所有数据集 - 支持分片文件")
    print("=" * 70)

    # 加载参考数据集
    print("步骤1: 加载 FlowSorted.Blood.450k 参考数据集")
    reference_beta, reference_pheno = load_flowsorted_reference()

    if reference_beta is None:
        print("❌ 无法加载参考数据集，将只进行健康样本筛选")
        # 创建空的参考数据，这样代码还能运行但跳过细胞比例计算
        reference_beta = None
        reference_pheno = None

    # 查找所有存在的文件
    available_datasets = []
    for dataset in target_datasets:
        pheno_file = os.path.join(data_dir, f"{dataset}_pheno.csv")

        # 检查是否存在任何形式的beta文件
        has_beta_files = False
        if os.path.exists(os.path.join(data_dir, f"{dataset}_beta.csv")):
            has_beta_files = True
        else:
            for i in range(1, 20):
                if os.path.exists(os.path.join(data_dir, f"{dataset}_beta_{i}.csv")):
                    has_beta_files = True
                    break

        if os.path.exists(pheno_file) and has_beta_files:
            available_datasets.append(dataset)
            print(f"✓ 找到完整文件: {dataset}")
        else:
            print(f"✗ 文件不完整: {dataset}")

    print(f"\n步骤2: 找到 {len(available_datasets)} 个完整的数据集")

    # 存储处理结果
    processing_results = {}
    cell_types_available = []

    # 处理每个数据集
    for dataset in available_datasets:
        result = process_single_dataset_with_cell_proportions(
            dataset, data_dir, pheno_output_dir, reference_beta, reference_pheno
        )
        processing_results[dataset] = result

    # 生成最终报告
    print(f"\n步骤5: 生成处理报告...")
    generate_comprehensive_report(processing_results, other_files_dir, reference_beta)

    return processing_results, cell_types_available, available_datasets, target_datasets, other_files_dir


def generate_comprehensive_report(processing_results, output_dir, reference_beta):
    """生成综合处理报告"""
    successful = [ds for ds, info in processing_results.items() if info.get('status') == 'success']
    failed = [ds for ds, info in processing_results.items() if info.get('status') == 'error']

    successful_with_cells = [ds for ds in successful if processing_results[ds].get('has_cell_proportions', False)]
    successful_without_cells = [ds for ds in successful if
                                not processing_results[ds].get('has_cell_proportions', False)]

    total_samples = sum(
        [info.get('sample_count', 0) for ds, info in processing_results.items() if info.get('status') == 'success'])
    total_original_samples = sum(
        [info.get('health_stats', {}).get('total_samples', 0) for ds, info in processing_results.items() if
         info.get('status') == 'success'])

    # 收集所有细胞类型
    all_cell_types = set()
    for ds in successful_with_cells:
        cell_types = processing_results[ds].get('cell_types', [])
        if cell_types:
            all_cell_types.update(cell_types)

    report = {
        'processing_summary': {
            'total_target_datasets': 17,
            'successful_datasets': len(successful),
            'failed_datasets': len(failed),
            'datasets_with_cell_proportions': len(successful_with_cells),
            'datasets_without_cell_proportions': len(successful_without_cells),
            'total_healthy_samples': total_samples,
            'total_original_samples': total_original_samples,
            'healthy_sample_ratio': round(total_samples / total_original_samples,
                                          3) if total_original_samples > 0 else 0,
            'all_cell_types': list(all_cell_types),
            'reference_data_available': reference_beta is not None
        },
        'reference_info': {
            'name': 'FlowSorted.Blood.450k',
            'method': 'Linear Regression Deconvolution',
            'file_path': '/mnt/sfs_turbo/wanghui/protein/file/FlowSorted.Blood.450k.csv'
        },
        'detailed_results': {
            ds: {
                'sample_count': processing_results[ds].get('sample_count', 0),
                'has_cell_proportions': processing_results[ds].get('has_cell_proportions', False),
                'cell_types': processing_results[ds].get('cell_types', []),
                'health_stats': processing_results[ds].get('health_stats', {}),
                'output_file': processing_results[ds].get('pheno_file', '')
            } for ds in successful
        },
        'failed_datasets': {
            ds: processing_results[ds].get('error_message', 'Unknown error') for ds in failed
        },
        'file_locations': {
            'pheno_files': '/mnt/sfs_turbo/wanghui/protein/data/GSE-data/',
            'report_files': output_dir
        }
    }

    report_path = os.path.join(output_dir, 'complete_cell_proportion_analysis_report.json')
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"✓ 综合报告已保存: {report_path}")

    # 打印最终摘要
    print(f"\n" + "=" * 70)
    print("处理完成摘要")
    print("=" * 70)
    print(f"目标数据集: 17 个")
    print(f"成功处理: {len(successful)} 个数据集")
    print(f"  - 包含细胞比例: {len(successful_with_cells)} 个")
    print(f"  - 无细胞比例: {len(successful_without_cells)} 个")
    print(f"处理失败: {len(failed)} 个数据集")
    print(f"总健康样本: {total_samples} 个 (从 {total_original_samples} 个原始样本中筛选)")
    print(f"健康样本比例: {report['processing_summary']['healthy_sample_ratio'] * 100:.1f}%")
    print(f"参考数据可用: {reference_beta is not None}")
    if all_cell_types:
        print(f"发现的细胞类型: {list(all_cell_types)}")
    print(f"\n文件位置:")
    print(f"- 表型文件: {report['file_locations']['pheno_files']}")
    print(f"- 报告文件: {report_path}")


# 执行完整处理流程
print("开始完整的细胞比例计算流程（支持分片文件）...")
try:
    processing_results, cell_types, available_datasets, target_datasets, other_files_dir = process_all_datasets_with_cell_proportions()
    print(f"\n🎉 所有处理完成！")

    # 打印成功的数据集
    successful = [ds for ds, info in processing_results.items() if info.get('status') == 'success']
    print(f"\n成功处理的数据集 ({len(successful)} 个): {successful}")

    if successful:
        print(f"\n每个数据集的样本数量:")
        for dataset in successful:
            info = processing_results[dataset]
            cell_status = "有细胞比例" if info.get('has_cell_proportions', False) else "无细胞比例"
            print(f"  - {dataset}: {info.get('sample_count', 0)} 个健康样本 ({cell_status})")

except Exception as e:
    print(f"\n❌ 处理过程中出现错误: {e}")
    import traceback

    traceback.print_exc()