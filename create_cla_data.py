import pandas as pd
import os
import re
from rdkit import Chem
from rdkit import rdBase
from functools import reduce
from sklearn.model_selection import train_test_split
import numpy as np

# 禁用RDKit的详细错误日志，以便我们自己捕获和处理
rdBase.DisableLog('rdApp.error')

# --- 核心辅助函数 ---

def diagnose_and_fix_smiles(smiles: str) -> dict:
    """
    诊断并尝试修复单个SMILES字符串。
    """
    if not isinstance(smiles, str) or not smiles:
        return {'status': 'unfixable', 'smiles': None, 'error': 'Input is not a valid string.'}

    # 1. 尝试原始SMILES
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is not None:
            return {'status': 'valid', 'smiles': smiles, 'error': None}
    except Exception as e:
        pass # 继续尝试修复

    # 2. 应用修复策略
    smiles_fixed = smiles
    
    # 处理常见的元素符号错误
    element_typo_replacements = {'IN': '[In]'}
    for typo, correction in element_typo_replacements.items():
        smiles_fixed = re.sub(r'(^|[^a-zA-Z])' + re.escape(typo) + r'($|[^a-zA-Z])', r'\1' + correction + r'\2', smiles_fixed)

    # 转换 [N+H] 风格为 [NH+] 风格
    charge_hydrogen_format_replacements = {
        '[N+H]': '[NH+]', '[N+H2]': '[NH2+]', '[N+H3]': '[NH3+]', '[n+H]': '[nH+]', '[O+H]': '[OH+]'
    }
    for old, new in charge_hydrogen_format_replacements.items():
        smiles_fixed = smiles_fixed.replace(old, new)
    
    # 移除其他带电原子中不必要的显式氢原子
    smiles_fixed = re.sub(r'\[([A-Za-z@\*]+[+\-]\d+)H\d*\]', r'[\1]', smiles_fixed)

    # 3. 尝试解析修复后的SMILES
    try:
        mol = Chem.MolFromSmiles(smiles_fixed, sanitize=True)
        if mol is not None:
            if smiles_fixed != smiles:
                return {'status': 'fixed', 'smiles': smiles_fixed, 'error': None}
            else:
                return {'status': 'valid', 'smiles': smiles, 'error': None}
    except Exception as e:
        return {'status': 'unfixable', 'smiles': None, 'error': str(e)}

    return {'status': 'unfixable', 'smiles': None, 'error': 'Initial parsing failed and no fixes were applicable.'}


def create_stratification_label(df, task_columns):
    """
    为分层采样创建分层标签
    对于多任务分类，创建一个组合标签用于分层
    """
    # 方法1：使用主要任务（有数据最多的任务）进行分层
    task_counts = {}
    for task in task_columns:
        task_counts[task] = df[task].notna().sum()
    
    # 选择数据最多的任务作为主分层任务
    main_task = max(task_counts, key=task_counts.get)
    print(f"  使用 {main_task} 作为主分层任务 ({task_counts[main_task]} 个有效样本)")
    
    # 创建分层标签
    stratify_labels = []
    for idx, row in df.iterrows():
        if pd.notna(row[main_task]):
            # 使用主任务的标签
            stratify_labels.append(f"{main_task}_{int(row[main_task])}")
        else:
            # 如果主任务没有标签，查找其他有标签的任务
            found_label = False
            for task in task_columns:
                if pd.notna(row[task]):
                    stratify_labels.append(f"{task}_{int(row[task])}")
                    found_label = True
                    break
            if not found_label:
                stratify_labels.append("no_label")
    
    return stratify_labels


def process_and_split_datasets(input_dir, output_dir, files_to_skip=None, split_ratios=(0.8, 0.1, 0.1), random_state=42):
    """
    加载、修复、合并并分割CSV文件，输出train.csv, val.csv, test.csv
    专门针对分类任务，使用分层采样
    """
    if files_to_skip is None:
        files_to_skip = []
    
    skip_list_lower = [item.lower().strip() for item in files_to_skip]
    all_processed_dfs = []

    print(f"开始处理目录 '{input_dir}' ...")
    print(f"分割比例: 训练集 {split_ratios[0]:.1%}, 验证集 {split_ratios[1]:.1%}, 测试集 {split_ratios[2]:.1%}")
    print(f"使用分层采样进行分类任务分割")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理所有任务文件
    for filename in os.listdir(input_dir):
        if not filename.endswith('.csv'): 
            continue
        
        task_name_original = filename.split('.')[0]
        task_name_normalized = task_name_original.strip().lower()

        if task_name_normalized in skip_list_lower: 
            print(f"跳过文件: {filename}")
            continue
            
        print(f"处理: {filename}")
        file_path = os.path.join(input_dir, filename)

        try:
            df = pd.read_csv(file_path)
            if 'SMILES' in df.columns: 
                df.rename(columns={'SMILES': 'smiles'}, inplace=True)
            if 'smiles' not in df.columns or 'Label' not in df.columns: 
                continue

            # 诊断和修复SMILES
            diagnostics = df['smiles'].apply(diagnose_and_fix_smiles)
            df['smiles_fixed'] = diagnostics.apply(lambda x: x['smiles'])
            df['status'] = diagnostics.apply(lambda x: x['status'])

            # 清洗数据：只保留有效和已修复的SMILES
            cleaned_df = df[df['status'].isin(['valid', 'fixed'])].copy()
            if cleaned_df.empty: 
                continue

            # 重命名Label列并创建最终数据框（不包含Inchikey）
            cleaned_df.rename(columns={'Label': task_name_original}, inplace=True)
            final_df = cleaned_df[['smiles_fixed', task_name_original]].copy()
            final_df.rename(columns={'smiles_fixed': 'smiles'}, inplace=True)
            
            if not final_df.empty:
                all_processed_dfs.append(final_df)
                print(f"  -> {len(final_df)} 个有效样本")
                
                # 显示分类标签分布
                label_counts = final_df[task_name_original].value_counts().sort_index()
                print(f"  -> 标签分布: {dict(label_counts)}")

        except Exception as e:
            print(f"  -> 处理文件 {filename} 时出错: {e}")

    # 合并所有任务数据
    if not all_processed_dfs:
        print("没有有效数据可以合并！")
        return

    print(f"\n合并 {len(all_processed_dfs)} 个任务的数据...")
    
    # 使用smiles作为唯一键进行合并（去掉Inchikey）
    complete_dataset = all_processed_dfs[0]
    for df in all_processed_dfs[1:]:
        complete_dataset = pd.merge(complete_dataset, df, on='smiles', how='outer')
    
    print(f"合并后数据集: {len(complete_dataset)} 个唯一分子")
    
    # 获取任务列（除了smiles列）
    task_columns = [col for col in complete_dataset.columns if col != 'smiles']
    print(f"任务列: {task_columns}")
    
    # 显示各任务的标签分布
    print(f"\n各任务标签分布:")
    for task in task_columns:
        valid_data = complete_dataset[task].dropna()
        if len(valid_data) > 0:
            label_counts = valid_data.value_counts().sort_index()
            print(f"  {task}: {dict(label_counts)} (总计 {len(valid_data)} 样本)")
        else:
            print(f"  {task}: 无有效数据")
    
    # 分割数据集（使用分层采样）
    total_samples = len(complete_dataset)
    test_size = split_ratios[2]
    val_size = split_ratios[1]
    
    print(f"\n开始分层分割数据集...")
    
    if total_samples < 10:
        print("警告: 数据量太少，不进行分割")
        train_data = complete_dataset.copy()
        val_data = pd.DataFrame(columns=complete_dataset.columns)
        test_data = pd.DataFrame(columns=complete_dataset.columns)
    else:
        try:
            # 创建分层标签
            stratify_labels = create_stratification_label(complete_dataset, task_columns)
            
            # 检查每个分层标签的样本数量
            from collections import Counter
            label_counts = Counter(stratify_labels)
            print(f"  分层标签分布: {dict(label_counts)}")
            
            # 过滤掉样本数量过少的标签（至少需要2个样本才能分层）
            valid_indices = []
            valid_stratify_labels = []
            
            for i, label in enumerate(stratify_labels):
                if label_counts[label] >= 2:  # 至少2个样本才能分层
                    valid_indices.append(i)
                    valid_stratify_labels.append(label)
            
            if len(valid_indices) < len(complete_dataset) * 0.5:
                print("  警告: 可用于分层的样本太少，改用随机分割")
                # 使用随机分割
                train_val_data, test_data = train_test_split(
                    complete_dataset, 
                    test_size=test_size, 
                    random_state=random_state
                )
                
                if len(train_val_data) < 5:
                    train_data = train_val_data.copy()
                    val_data = pd.DataFrame(columns=complete_dataset.columns)
                else:
                    val_size_adjusted = val_size / (1 - test_size)
                    train_data, val_data = train_test_split(
                        train_val_data,
                        test_size=val_size_adjusted,
                        random_state=random_state
                    )
            else:
                # 使用分层分割
                valid_data = complete_dataset.iloc[valid_indices]
                remaining_data = complete_dataset.drop(complete_dataset.index[valid_indices])
                
                print(f"  使用分层分割: {len(valid_data)} 样本")
                print(f"  无法分层的样本: {len(remaining_data)} 样本")
                
                # 对可分层的数据进行分层分割
                train_val_valid, test_valid = train_test_split(
                    valid_data,
                    test_size=test_size,
                    stratify=valid_stratify_labels,
                    random_state=random_state
                )
                
                # 对剩余数据进行随机分割
                if len(remaining_data) > 0:
                    train_val_remain, test_remain = train_test_split(
                        remaining_data,
                        test_size=test_size,
                        random_state=random_state
                    )
                    
                    # 合并分层和随机分割的结果
                    train_val_data = pd.concat([train_val_valid, train_val_remain], ignore_index=True)
                    test_data = pd.concat([test_valid, test_remain], ignore_index=True)
                else:
                    train_val_data = train_val_valid
                    test_data = test_valid
                
                # 继续分割训练集和验证集
                if len(train_val_data) < 5:
                    train_data = train_val_data.copy()
                    val_data = pd.DataFrame(columns=complete_dataset.columns)
                else:
                    val_size_adjusted = val_size / (1 - test_size)
                    
                    # 尝试对训练+验证集进行分层分割
                    try:
                        train_val_stratify = create_stratification_label(train_val_data, task_columns)
                        train_val_label_counts = Counter(train_val_stratify)
                        
                        # 检查是否可以分层
                        can_stratify = all(count >= 2 for count in train_val_label_counts.values())
                        
                        if can_stratify and len(set(train_val_stratify)) > 1:
                            train_data, val_data = train_test_split(
                                train_val_data,
                                test_size=val_size_adjusted,
                                stratify=train_val_stratify,
                                random_state=random_state
                            )
                            print("  训练/验证集也使用了分层分割")
                        else:
                            train_data, val_data = train_test_split(
                                train_val_data,
                                test_size=val_size_adjusted,
                                random_state=random_state
                            )
                            print("  训练/验证集使用随机分割")
                    except:
                        train_data, val_data = train_test_split(
                            train_val_data,
                            test_size=val_size_adjusted,
                            random_state=random_state
                        )
                        print("  训练/验证集分层失败，使用随机分割")
            
        except Exception as e:
            print(f"  分层分割失败: {e}")
            print("  改用随机分割...")
            
            # 降级到随机分割
            train_val_data, test_data = train_test_split(
                complete_dataset, 
                test_size=test_size, 
                random_state=random_state
            )
            
            if len(train_val_data) < 5:
                train_data = train_val_data.copy()
                val_data = pd.DataFrame(columns=complete_dataset.columns)
            else:
                val_size_adjusted = val_size / (1 - test_size)
                train_data, val_data = train_test_split(
                    train_val_data,
                    test_size=val_size_adjusted,
                    random_state=random_state
                )

    print(f"\n分割结果:")
    print(f"  训练集: {len(train_data):,} 样本")
    print(f"  验证集: {len(val_data):,} 样本")
    print(f"  测试集: {len(test_data):,} 样本")

    # 显示各数据集的标签分布
    for dataset_name, dataset in [("训练集", train_data), ("验证集", val_data), ("测试集", test_data)]:
        if not dataset.empty:
            print(f"\n{dataset_name}标签分布:")
            for task in task_columns:
                valid_data = dataset[task].dropna()
                if len(valid_data) > 0:
                    label_counts = valid_data.value_counts().sort_index()
                    print(f"  {task}: {dict(label_counts)}")

    # 保存分割后的数据集
    if not train_data.empty:
        train_path = os.path.join(output_dir, 'train.csv')
        train_data.to_csv(train_path, index=False)
        print(f"\n训练集已保存至: {train_path}")

    if not val_data.empty:
        val_path = os.path.join(output_dir, 'val.csv')
        val_data.to_csv(val_path, index=False)
        print(f"验证集已保存至: {val_path}")

    if not test_data.empty:
        test_path = os.path.join(output_dir, 'test.csv')
        test_data.to_csv(test_path, index=False)
        print(f"测试集已保存至: {test_path}")

    print(f"\n处理完成！输出文件:")
    print(f"  📁 {output_dir}/train.csv")
    print(f"  📁 {output_dir}/val.csv")
    print(f"  📁 {output_dir}/test.csv")


# --- 执行 ---
if __name__ == '__main__':
    # 测试诊断功能
    print("--- 测试诊断功能 ---")
    test_smiles = 'BrN([C-2]C)([C-2]C)([C-2]C)[C-2]C'
    result = diagnose_and_fix_smiles(test_smiles)
    print(f"原始SMILES: {test_smiles}")
    print(f"诊断结果: {result}")
    print("-" * 50)
    
    # 运行数据分割和合并流程
    INPUT_CLA_DIR = 'ademt_data/wash_cla'  # 分类数据目录
    OUTPUT_DIR = 'ademt_data/cla_3'
    SKIP_FILES = []  # 可以指定要跳过的文件
    
    # 设置分割参数
    SPLIT_RATIOS = (0.8, 0.1, 0.1)  # 训练集, 验证集, 测试集
    RANDOM_STATE = 587 # 随机种子

    process_and_split_datasets(
        input_dir=INPUT_CLA_DIR, 
        output_dir=OUTPUT_DIR, 
        files_to_skip=SKIP_FILES,
        split_ratios=SPLIT_RATIOS,
        random_state=RANDOM_STATE
    )     # 训练集, 验证集, 测试集

    process_and_split_datasets(
        input_dir=INPUT_CLA_DIR, 
        output_dir=OUTPUT_DIR, 
        files_to_skip=SKIP_FILES,
        split_ratios=SPLIT_RATIOS,
        random_state=RANDOM_STATE
    )