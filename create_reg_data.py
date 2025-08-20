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

# --- 核心辅助函数 (带有诊断功能) ---

def smiles_to_inchikey(smiles_list):
    """将SMILES列表转换为InChIKey列表。"""
    inchikeys = []
    for smiles in smiles_list:
        if pd.isna(smiles):
            inchikeys.append(None)
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            inchikey = Chem.MolToInchiKey(mol)
            inchikeys.append(inchikey)
        else:
            inchikeys.append(None)
    return inchikeys


def diagnose_and_fix_smiles(smiles: str) -> dict:
    """
    诊断并尝试修复单个SMILES字符串。

    Returns:
        一个包含诊断结果的字典:
        {'status': 'valid'/'fixed'/'unfixable', 'smiles': str/None, 'error': str/None}
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

    # 2. 应用一系列修复策略
    smiles_fixed = smiles
    
    # 策略 A: 处理常见的元素符号错误
    element_typo_replacements = {'IN': '[In]'}
    for typo, correction in element_typo_replacements.items():
        smiles_fixed = re.sub(r'(^|[^a-zA-Z])' + re.escape(typo) + r'($|[^a-zA-Z])', r'\1' + correction + r'\2', smiles_fixed)

    # 策略 B: 转换 [N+H] 风格为 [NH+] 风格
    charge_hydrogen_format_replacements = {
        '[N+H]': '[NH+]', '[N+H2]': '[NH2+]', '[N+H3]': '[NH3+]', '[n+H]': '[nH+]', '[O+H]': '[OH+]'
    }
    for old, new in charge_hydrogen_format_replacements.items():
        smiles_fixed = smiles_fixed.replace(old, new)
    
    # 策略 C: 移除其他带电原子中不必要的显式氢原子
    smiles_fixed = re.sub(r'\[([A-Za-z@\*]+[+\-]\d+)H\d*\]', r'[\1]', smiles_fixed)

    # 3. 尝试解析修复后的SMILES并捕获错误
    try:
        mol = Chem.MolFromSmiles(smiles_fixed, sanitize=True)
        if mol is not None:
            # 检查修复是否真的改变了字符串
            if smiles_fixed != smiles:
                return {'status': 'fixed', 'smiles': smiles_fixed, 'error': None}
            else:
                # 这种情况理论上不应该发生，因为原始的已经检查过了
                return {'status': 'valid', 'smiles': smiles, 'error': None}
    except Exception as e:
        # 如果修复后仍然失败，捕获错误信息
        return {'status': 'unfixable', 'smiles': None, 'error': str(e)}

    # 如果修复后没有变化且原始SMILES无效
    return {'status': 'unfixable', 'smiles': None, 'error': 'Initial parsing failed and no fixes were applicable.'}


def split_dataset(df, task_name, test_size=0.1, val_size=0.1, random_state=42):
    """
    将数据集按照8:1:1的比例分割为训练集、验证集和测试集
    
    Args:
        df: 包含数据的DataFrame
        task_name: 任务名称
        test_size: 测试集比例
        val_size: 验证集比例
        random_state: 随机种子
    
    Returns:
        train_df, val_df, test_df: 分割后的三个数据集
    """
    if len(df) < 10:  # 如果数据太少，不进行分割
        print(f"  警告: {task_name} 数据量太少 ({len(df)} 样本)，不进行分割")
        return df.copy(), pd.DataFrame(), pd.DataFrame()
    
    # 首先分出测试集
    train_val_df, test_df = train_test_split(
        df, 
        test_size=test_size, 
        random_state=random_state,
        stratify=None  # 对于回归任务不使用分层抽样
    )
    
    # 再从训练+验证集中分出验证集
    if len(train_val_df) < 5:  # 如果剩余数据太少
        return train_val_df.copy(), pd.DataFrame(), test_df
    
    # 计算验证集在剩余数据中的比例
    val_size_adjusted = val_size / (1 - test_size)
    
    train_df, val_df = train_test_split(
        train_val_df,
        test_size=val_size_adjusted,
        random_state=random_state
    )
    
    return train_df, val_df, test_df


def process_and_split_datasets(input_dir, output_dir, files_to_skip=None, split_ratios=(0.8, 0.1, 0.1), random_state=42, report_unfixable=True):
    """
    加载、修复、诊断、清洗、分割并合并CSV文件。
    
    Args:
        input_dir: 输入目录路径
        output_dir: 输出目录路径
        files_to_skip: 要跳过的文件列表
        split_ratios: 分割比例 (train, val, test)
        random_state: 随机种子
        report_unfixable: 是否报告无法修复的SMILES
    """
    if files_to_skip is None:
        files_to_skip = []
    
    skip_list_lower = [item.lower().strip() for item in files_to_skip]
    
    # 存储所有任务的分割数据
    all_train_dfs = []
    all_val_dfs = []
    all_test_dfs = []
    
    unfixable_report = []
    split_summary = []

    print(f"开始处理目录 '{input_dir}' ...")
    print(f"分割比例: 训练集 {split_ratios[0]:.1%}, 验证集 {split_ratios[1]:.1%}, 测试集 {split_ratios[2]:.1%}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in os.listdir(input_dir):
        if not filename.endswith('.csv'): 
            continue
        
        task_name_original = filename.split('.')[0]
        task_name_normalized = task_name_original.strip().lower()

        if task_name_normalized in skip_list_lower: 
            print(f"跳过文件: {filename}")
            continue
            
        print(f"\n--- 正在处理: {filename} ---")
        file_path = os.path.join(input_dir, filename)

        try:
            df = pd.read_csv(file_path)
            if 'SMILES' in df.columns: 
                df.rename(columns={'SMILES': 'smiles'}, inplace=True)
            if 'smiles' not in df.columns or 'Label' not in df.columns: 
                print(f"  跳过: 缺少必要的列 (smiles 或 Label)")
                continue

            print(f"  原始数据: {len(df)} 样本")

            # 应用诊断和修复函数
            diagnostics = df['smiles'].apply(diagnose_and_fix_smiles)
            df['smiles_fixed'] = diagnostics.apply(lambda x: x['smiles'])
            df['status'] = diagnostics.apply(lambda x: x['status'])
            df['error_msg'] = diagnostics.apply(lambda x: x['error'])

            # 报告无法修复的SMILES
            if report_unfixable:
                unfixable_df = df[df['status'] == 'unfixable']
                for _, row in unfixable_df.iterrows():
                    unfixable_report.append({
                        'file': filename,
                        'original_smiles': row['smiles'],
                        'error': row['error_msg']
                    })

            # 清洗数据
            cleaned_df = df.dropna(subset=['smiles_fixed']).copy()
            if cleaned_df.empty: 
                print(f"  清洗后无有效数据")
                continue

            print(f"  清洗后数据: {len(cleaned_df)} 样本")

            # 生成InChIKey和重命名
            cleaned_df['Inchikey'] = smiles_to_inchikey(cleaned_df['smiles_fixed'])
            cleaned_df.rename(columns={'Label': task_name_original}, inplace=True)
            final_df = cleaned_df[['smiles_fixed', 'Inchikey', task_name_original]].copy()
            final_df.rename(columns={'smiles_fixed': 'smiles'}, inplace=True)
            
            # 移除InChIKey生成失败的行
            final_df = final_df.dropna(subset=['Inchikey'])
            
            if final_df.empty:
                print(f"  生成InChIKey后无有效数据")
                continue

            print(f"  最终数据: {len(final_df)} 样本")

            # 分割数据集
            train_df, val_df, test_df = split_dataset(
                final_df, 
                task_name_original, 
                test_size=split_ratios[2], 
                val_size=split_ratios[1], 
                random_state=random_state
            )

            print(f"  分割结果: 训练集 {len(train_df)}, 验证集 {len(val_df)}, 测试集 {len(test_df)}")

            # 记录分割统计
            split_summary.append({
                'task': task_name_original,
                'total': len(final_df),
                'train': len(train_df),
                'val': len(val_df),
                'test': len(test_df),
                'train_pct': len(train_df) / len(final_df) * 100 if len(final_df) > 0 else 0,
                'val_pct': len(val_df) / len(final_df) * 100 if len(final_df) > 0 else 0,
                'test_pct': len(test_df) / len(final_df) * 100 if len(final_df) > 0 else 0
            })

            # 保存单个任务的分割结果
            task_output_dir = os.path.join(output_dir, 'individual_tasks', task_name_original)
            os.makedirs(task_output_dir, exist_ok=True)
            
            if not train_df.empty:
                train_df.to_csv(os.path.join(task_output_dir, 'train.csv'), index=False)
            if not val_df.empty:
                val_df.to_csv(os.path.join(task_output_dir, 'val.csv'), index=False)
            if not test_df.empty:
                test_df.to_csv(os.path.join(task_output_dir, 'test.csv'), index=False)

            # 添加到总的数据集列表
            if not train_df.empty:
                all_train_dfs.append(train_df)
            if not val_df.empty:
                all_val_dfs.append(val_df)
            if not test_df.empty:
                all_test_dfs.append(test_df)

        except Exception as e:
            print(f"  -> 处理文件 {filename} 时发生严重错误: {e}")

    # 合并所有任务的数据集
    print(f"\n--- 合并所有任务的数据集 ---")
    
    def merge_datasets(df_list, set_name):
        """合并数据集列表"""
        if not df_list:
            print(f"  {set_name}: 没有数据可合并")
            return pd.DataFrame()
        
        print(f"  {set_name}: 合并 {len(df_list)} 个任务的数据")
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, on=['smiles', 'Inchikey'], how='outer'), 
            df_list
        )
        print(f"    合并后: {len(merged_df)} 个唯一分子")
        return merged_df

    # 合并训练集
    merged_train = merge_datasets(all_train_dfs, "训练集")
    if not merged_train.empty:
        train_output_path = os.path.join(output_dir, 'train.csv')
        merged_train.drop(['Inchikey'], axis=1, inplace=True, errors='ignore')  # 移除InChIKey列
        merged_train.to_csv(train_output_path, index=False)
        print(f"    保存至: {train_output_path}")

    # 合并验证集
    merged_val = merge_datasets(all_val_dfs, "验证集")
    if not merged_val.empty:
        val_output_path = os.path.join(output_dir, 'val.csv')
        merged_val.drop(['Inchikey'], axis=1, inplace=True, errors='ignore')  # 移除InChIKey列
        merged_val.to_csv(val_output_path, index=False)
        print(f"    保存至: {val_output_path}")

    # 合并测试集
    merged_test = merge_datasets(all_test_dfs, "测试集")
    if not merged_test.empty:
        test_output_path = os.path.join(output_dir, 'test.csv')
        merged_test.drop(['Inchikey'], axis=1, inplace=True, errors='ignore')  # 移除InChIKey列
        merged_test.to_csv(test_output_path, index=False)
        print(f"    保存至: {test_output_path}")

    # 保存分割统计报告
    if split_summary:
        summary_df = pd.DataFrame(split_summary)
        summary_path = os.path.join(output_dir, 'split_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        
        print(f"\n--- 分割统计报告 ---")
        print(summary_df.to_string(index=False))
        print(f"详细统计已保存至: {summary_path}")

    # 合并所有数据（不分割）用于对比
    print(f"\n--- 生成完整合并数据集 ---")
    if all_train_dfs or all_val_dfs or all_test_dfs:
        # 收集所有任务的完整数据
        all_complete_dfs = []
        
        for filename in os.listdir(input_dir):
            if not filename.endswith('.csv'): 
                continue
            
            task_name_original = filename.split('.')[0]
            task_name_normalized = task_name_original.strip().lower()

            if task_name_normalized in skip_list_lower: 
                continue
            
            # 重新读取并处理（这次不分割）
            file_path = os.path.join(input_dir, filename)
            try:
                df = pd.read_csv(file_path)
                if 'SMILES' in df.columns: 
                    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
                if 'smiles' not in df.columns or 'Label' not in df.columns: 
                    continue

                # 诊断和修复
                diagnostics = df['smiles'].apply(diagnose_and_fix_smiles)
                df['smiles_fixed'] = diagnostics.apply(lambda x: x['smiles'])
                cleaned_df = df.dropna(subset=['smiles_fixed']).copy()
                
                if cleaned_df.empty:
                    continue

                # 处理最终数据
                cleaned_df['Inchikey'] = smiles_to_inchikey(cleaned_df['smiles_fixed'])
                cleaned_df.rename(columns={'Label': task_name_original}, inplace=True)
                final_df = cleaned_df[['smiles_fixed', 'Inchikey', task_name_original]].copy()
                final_df.rename(columns={'smiles_fixed': 'smiles'}, inplace=True)
                final_df = final_df.dropna(subset=['Inchikey'])
                
                if not final_df.empty:
                    all_complete_dfs.append(final_df)
                    
            except Exception as e:
                continue

        if all_complete_dfs:
            complete_merged = reduce(
                lambda left, right: pd.merge(left, right, on=['smiles', 'Inchikey'], how='outer'), 
                all_complete_dfs
            )
            complete_path = os.path.join(output_dir, 'complete_dataset.csv')
            complete_merged.to_csv(complete_path, index=False)
            print(f"完整数据集: {len(complete_merged)} 个分子，保存至: {complete_path}")

    # 打印无法修复的SMILES报告
    if unfixable_report:
        print("\n--- 无法修复的SMILES报告 ---")
        report_df = pd.DataFrame(unfixable_report)
        unfixable_path = os.path.join(output_dir, 'unfixable_smiles_report.csv')
        report_df.to_csv(unfixable_path, index=False)
        print(f"发现 {len(unfixable_report)} 个无法修复的SMILES")
        print(f"详细报告已保存至: {unfixable_path}")
        
        if len(unfixable_report) <= 10:
            print(report_df.to_string(index=False))

    print(f"\n=== 处理完成 ===")
    print(f"所有结果保存在目录: {output_dir}")
    print(f"  - train.csv: 合并的训练集")
    print(f"  - val.csv: 合并的验证集") 
    print(f"  - test.csv: 合并的测试集")
    print(f"  - complete_dataset.csv: 完整合并数据集")
    print(f"  - individual_tasks/: 各任务的单独分割结果")
    print(f"  - split_summary.csv: 分割统计报告")


# --- 执行 ---
if __name__ == '__main__':
    # --- 测试新的诊断函数 ---
    print("--- 测试诊断功能 ---")
    test_smiles = 'BrN([C-2]C)([C-2]C)([C-2]C)[C-2]C'
    result = diagnose_and_fix_smiles(test_smiles)
    print(f"原始SMILES: {test_smiles}")
    print(f"诊断结果: {result}")
    print("-" * 50)
    
    # --- 运行数据分割和合并流程 ---
    INPUT_REG_DIR = 'ademt_data/wash_reg'
    OUTPUT_DIR = 'ademt_data/reg_3'
    SKIP_FILES = []  # 可以指定要跳过的文件，例如 ['logD', 'logP']
    
    # 设置分割参数
    SPLIT_RATIOS = (0.8, 0.1, 0.1)  # 训练集, 验证集, 测试集
    RANDOM_STATE = 50  # 随机种子，确保结果可重现

    process_and_split_datasets(
        input_dir=INPUT_REG_DIR, 
        output_dir=OUTPUT_DIR, 
        files_to_skip=SKIP_FILES,
        split_ratios=SPLIT_RATIOS,
        random_state=RANDOM_STATE,
        report_unfixable=True
    )