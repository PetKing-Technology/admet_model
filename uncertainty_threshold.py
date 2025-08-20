import pandas as pd
import numpy as np
import chemprop
import os
import json
from tqdm import tqdm

def get_predictions_and_uncertainty(smiles_list, model_path, preds_path='preds.csv'):
    """
    使用chemprop对SMILES列表进行预测，并返回预测值和不确定度。
    这是一个简化的函数，用于一次性获取所有任务的预测结果。

    参数:
    smiles_list (list): SMILES字符串列表。
    model_path (str): 训练好的chemprop模型文件路径 (.pt)。
    preds_path (str): 临时存储预测结果的路径。

    返回:
    tuple: (preds_df, uncertainty_df) 两个pandas DataFrame。
    """
    # 创建一个临时的SMILES文件，因为chemprop的输入需要一个文件路径
    temp_smiles_file = 'temp_smiles_for_pred.csv'
    pd.DataFrame({'smiles': smiles_list}).to_csv(temp_smiles_file, index=False)

    arguments = [
        '--test_path', temp_smiles_file,
        '--preds_path', preds_path,
        '--checkpoint_dir', model_path,
        '--num_workers', '0',
        '--uncertainty_method', 'ensemble',
    ]

    args = chemprop.args.PredictArgs().parse_args(arguments)
    preds, uncertainty = chemprop.train.make_predictions(
        args=args, 
        smiles=[[s] for s in smiles_list], 
        return_uncertainty=True
    )

    # 从chemprop的输出中读取列名
    pred_df_from_file = pd.read_csv(preds_path)
    task_names = [col for col in pred_df_from_file.columns if col != 'smiles']
    task_names = [
        'logp', 'pka_acidic', 'cl-plasma', 'logs', 'logvdss', 
        'pka_basic', 't12', 'logd', 'mdck', 'caco2', 'ppb', 'cl-int'
    ]
    preds_df = pd.DataFrame(preds, columns=task_names, index=smiles_list)
    uncertainty_df = pd.DataFrame(uncertainty, columns=task_names, index=smiles_list)
    
    return preds_df, uncertainty_df


def find_optimal_uncertainty_thresholds(
    test_set_path: str,
    model_path: str,
    tasks: list
) -> dict:
    """
    为多任务回归模型中的每个任务找到最优的不确定度阈值。

    参数:
    test_set_path (str): 测试集的CSV文件路径，第一列应为'smiles'，其余为任务的真实值。
    model_path (str): 训练好的chemprop多任务模型 (.pt) 的路径。
    tasks (list): 需要计算阈值的任务名称列表。

    返回:
    dict: 一个字典，键是任务名称，值是该任务的最佳不确定度阈值。
    """
    print("1. 加载测试数据...")
    test_df = pd.read_csv(test_set_path)
    smiles_list = test_df['smiles'].tolist()

    print("2. 使用Chemprop进行预测并计算不确定度...")
    # 假设你的多任务模型被保存在一个.pt文件里
    predictions_df, uncertainty_df = get_predictions_and_uncertainty(smiles_list, model_path)

    optimal_thresholds = {}

    print("3. 为每个任务计算最佳不确定度阈值...")
    for task in tqdm(tasks, desc="处理任务"):
        # --- 数据准备 ---
        # 合并真实值、预测值和不确定度
        task_data = pd.DataFrame({
            'true': test_df.set_index('smiles')[task],
            'pred': predictions_df[task],
            'uncertainty': uncertainty_df[task]
        }).dropna() # 非常重要：移除缺少真实值的样本

        if task_data.empty:
            print(f"任务 '{task}' 没有有效的测试数据，已跳过。")
            optimal_thresholds[task] = None
            continue

        # --- 定义“高/低准确度”的黄金标准 ---
        task_data['abs_error'] = (task_data['true'] - task_data['pred']).abs()
        error_threshold = task_data['abs_error'].median() # 使用中位绝对误差作为标准
        
        # is_low_accuracy = True (Positive), is_high_accuracy = False (Negative)
        task_data['is_low_accuracy'] = task_data['abs_error'] >= error_threshold

        # --- 遍历不确定度阈值，计算约登指数 ---
        # 使用排序后的唯一不确定度值作为候选阈值
        candidate_thresholds = sorted(task_data['uncertainty'].unique())
        best_j = -1
        best_threshold = -1

        for u_thresh in candidate_thresholds:
            # 预测为Positive (低准确度)
            predicted_positive = task_data['uncertainty'] > u_thresh
            # 预测为Negative (高准确度)
            predicted_negative = ~predicted_positive

            # 真实的Positive (的确是低准确度)
            actual_positive = task_data['is_low_accuracy']
            # 真实的Negative (的确是高准确度)
            actual_negative = ~actual_positive

            # 计算混淆矩阵
            TP = (predicted_positive & actual_positive).sum()
            TN = (predicted_negative & actual_negative).sum()
            FP = (predicted_positive & actual_negative).sum()
            FN = (predicted_negative & actual_positive).sum()

            # 计算敏感性和特异性
            sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
            specificity = TN / (TN + FP) if (TN + FP) > 0 else 0

            # 计算约登指数
            j_index = sensitivity + specificity - 1

            if j_index > best_j:
                best_j = j_index
                best_threshold = u_thresh
        
        optimal_thresholds[task] = best_threshold
        
    print("\n4. 计算完成！")
    return optimal_thresholds

# --- 如何使用 ---
if __name__ == '__main__':
    # 请替换为您的实际文件路径
    TEST_SET_FILE = 'ademt_data/reg_2/test.csv'  # 您的测试集文件
    MODEL_FILE = 'model_path/reg_model_nordkit2d_2/fold_0' # 您的多任务模型文件

    # 定义您的模型中的所有回归任务
    TASKS_LIST = [
        'logp', 'pka_acidic', 'cl-plasma', 'logs', 'logvdss', 
        'pka_basic', 't12', 'logd', 'mdck', 'caco2', 'ppb', 'cl-int'
    ]

    # 检查文件是否存在
    if not os.path.exists(TEST_SET_FILE) or not os.path.exists(MODEL_FILE):
        print("="*60)
        print("示例代码运行说明:")
        print(f"请将 'TEST_SET_FILE' 变量修改为您的测试集文件路径 ('{TEST_SET_FILE}')")
        print(f"请将 'MODEL_FILE' 变量修改为您的模型文件路径 ('{MODEL_FILE}')")
        print("这是一个示例，所以不会实际运行，除非您提供了正确的文件路径。")
        print("="*60)
        # 创建一个空字典作为示例输出
        final_thresholds = {task: "N/A" for task in TASKS_LIST}
    else:
        # 运行主函数
        final_thresholds = find_optimal_uncertainty_thresholds(
            test_set_path=TEST_SET_FILE,
            model_path=MODEL_FILE,
            tasks=TASKS_LIST
        )
    with open('uncertainty_thresholds.json', 'w') as f:
        json.dump(final_thresholds, f, indent=4)

    print("\n各终点对应的最优不确定度阈值:")
    print("-" * 40)
    for task, threshold in final_thresholds.items():
        if isinstance(threshold, float):
            print(f"{task:<12}: {threshold:.6f}")
        else:
            print(f"{task:<12}: {threshold}")
    print("-" * 40)