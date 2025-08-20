import os
import json
import time
import pandas as pd
import chemprop
import math # 引入math库以处理可能的无效数字
from scripts.Process_runner import  ProcessRunner
# --- 1. 配置区域 ---
MODEL_PATHS = {
    'reg': "/mnt/newdisk/fuli/ADMET/model_path/reg_model_nordkit2d_2/fold_0",
    'cla': "/mnt/newdisk/fuli/ADMET/model_path/cla_model_nordkit2d/fold_0"
}

UNCERTAINTY_THRESHOLDS = {
    "logp": 0.004031666415146518, "pka_acidic": 0.09403022325165011,
    "cl-plasma": 0.3946969731638035, "logs": 0.01319636964331572,
    "logvdss": 0.002377516110298794, "pka_basic": 0.04426391422279252,
    "t12": 0.017738629737266157, "logd": 0.010377547196687154,
    "mdck": 0.001470032041797964, "caco2": 0.003098651194555657,
    "ppb": 2.8860414398623107, "cl-int": 0.0039910197139036985
}
COLNAMES_DICT = {
   'cla': ['ames', 'lm-human', 'herg-10um', 'lm-mouse', 'f50', 'pgp_inh',
            'dili', 'herg', 'oatp1b1', 'cyp3a4-inh', 'cyp2d6-sub',
            'cyp3a4-sub', 'aggregators', 'oatp1b3', 'pgp_sub',
            'cyp2d6-inh', 'carcinogenicity', 'bbb'],
   'reg': ['logp', 'pka_acidic', 'cl-plasma', 'logs', 'logvdss',
            'pka_basic', 't12', 'logd', 'mdck', 'caco2', 'ppb', 'cl-int']
}

# --- 2. 辅助工具 & 函数 ---

class suppress_stdout_stderr(object):
    """用于屏蔽 chemprop 在终端的大量输出，使界面更整洁"""
    def __init__(self):
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        self.save_fds = (os.dup(1), os.dup(2))
    def __enter__(self):
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)
    def __exit__(self, *_):
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])

def classify_reg_confidence(uncertainty_value, threshold):
    """根据不确定度阈值判断回归任务的置信度"""
    if not isinstance(uncertainty_value, (float, int)):
        return "Invalid Input"
    if uncertainty_value > threshold:
        return "Low-confidence"
    else:
        return "High-confidence"

def classify_cla_confidence(probability):
    """根据预测概率判断分类任务的置信度"""
    if not isinstance(probability, (float, int)):
        return "Invalid Input"
    if probability < 0.3:
        return "excellent"
    elif probability < 0.7:
        return "medium"
    else:
        return "poor"

# *** 新增函数：用于解释回归预测值 ***
def get_regression_interpretation(task_name, value):
    """
    根据任务名称和预测值，返回一个解释性的字符串。
    """
    # 检查输入值是否有效
    if not isinstance(value, (int, float)) or math.isnan(value):
        return None

    interpretation = None
    try:
        if task_name == 'logp':
            interpretation = "excellent" if 0 <= value <= 3 else "poor"
        
        elif task_name == 'logs':
            interpretation = "excellent" if -4 <= value <= 0.5 else "poor"

        elif task_name == 'logd':  # 假设这是 logD at pH 7.4
            interpretation = "excellent" if 1 <= value <= 3 else "poor"

        elif task_name == 'logvdss':
            vdss = 10**value  # 从 logVdss 转换为 Vdss (L/kg)
            interpretation = "excellent" if 0.04 <= vdss <= 20 else "poor"
        
        elif task_name == 't12':
            if value > 8:
                interpretation = "excellent"
            elif 1 <= value <= 8:
                interpretation = "medium"
            else:  # value < 1
                interpretation = "poor"

        elif task_name == 'mdck':  # 假设模型输出单位为 10^-6 cm/s
            interpretation = "excellent" if value > 2 else "poor"
        
        elif task_name == 'caco2': # 假设模型输出为 logPapp
            interpretation = "excellent" if value > -5.15 else "poor"
        
        elif task_name == 'ppb':  # 假设模型输出为百分比 (%)
            interpretation = "excellent" if value <= 90 else "poor"
                
    except (TypeError, ValueError):
        # 如果计算出错（例如传入非数字），返回 None
        interpretation = None
        
    return interpretation

# --- 3. 核心预测函数 ---
# (此部分无变化)
def run_regression_predictions(smiles_list, smiles_file, pred_file):
    print("Running regression model predictions...")
    arguments = [
        '--test_path', smiles_file, '--preds_path', pred_file,
        '--checkpoint_dir', MODEL_PATHS['reg'], '--num_workers', '0',
        '--uncertainty_method', 'ensemble',
    ]
    args = chemprop.args.PredictArgs().parse_args(arguments)
    preds, un = chemprop.train.make_predictions(
        args=args, smiles=[[smi] for smi in smiles_list], return_uncertainty=True)
    task_names = COLNAMES_DICT['reg']
    preds_df = pd.DataFrame(preds, columns=task_names, index=smiles_list)
    un_df = pd.DataFrame(un, columns=[f"{col}_uncertainty" for col in task_names], index=smiles_list)
    return pd.concat([preds_df, un_df], axis=1)

def run_classification_predictions(smiles_list, smiles_file, pred_file):
    print("Running classification model predictions...")
    arguments = [
        '--test_path', smiles_file, '--preds_path', pred_file,
        '--checkpoint_dir', MODEL_PATHS['cla'], '--num_workers', '0',
    ]
    args = chemprop.args.PredictArgs().parse_args(arguments)
    preds = chemprop.train.make_predictions(args=args, smiles=[[smi] for smi in smiles_list])
    preds_df = pd.DataFrame(preds, columns=COLNAMES_DICT['cla'], index=smiles_list)
    return preds_df

# --- 4. 主流程函数 ---

def main(smiles_list):
    """
    主函数，协调回归和分类预测，并按要求格式化输出。
    """
    with ProcessRunner() as runner:
        smiles_file = 'temp_input.csv'
        pred_file = 'temp_output.csv'
        pd.DataFrame({'smiles': smiles_list}).to_csv(smiles_file, index=False)

        reg_results_df = run_regression_predictions(smiles_list, smiles_file, pred_file)
        cla_results_df = run_classification_predictions(smiles_list, smiles_file, pred_file)


    final_json_output = {}
    for smiles in smiles_list:
        final_json_output[smiles] = {}

        # *** 修改部分：处理回归结果并添加interpretation ***
        for task in COLNAMES_DICT['reg']:
            pred_val = reg_results_df.loc[smiles, task]
            un_val = reg_results_df.loc[smiles, f"{task}_uncertainty"]
            
            # 1. 创建基础的结果字典
            task_output = {
                "prediction": float(pred_val),
                "uncertainty": float(un_val),
                "confidence": classify_reg_confidence(un_val, UNCERTAINTY_THRESHOLDS.get(task))
            }

            # 2. 调用新函数获取解释，如果存在则添加到字典中
            interpretation = get_regression_interpretation(task, pred_val)
            if interpretation:
                task_output["interpretation"] = interpretation

            # 3. 将完整的任务结果存入最终输出
            final_json_output[smiles][task] = task_output

        # 处理分类结果 (无变化)
        for task in COLNAMES_DICT['cla']:
            prob_val = cla_results_df.loc[smiles, task]
            final_json_output[smiles][task] = {
                "prediction": float(prob_val),
                "confidence": classify_cla_confidence(prob_val)
            }
            
    return final_json_output

# --- 5. 脚本执行入口 ---
# (此部分无变化)
if __name__ == '__main__':
    test_smiles_list = ["CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"]
    print("Prediction process started...")
    start_time = time.time()
    with suppress_stdout_stderr():
        final_predictions = main(test_smiles_list)
    end_time = time.time()
    print(f"\nPrediction finished in {end_time - start_time:.2f} seconds.")

    print(json.dumps(final_predictions, indent=4))
    with open('prediction_output.json', 'w') as f:
        json.dump(final_predictions, f, indent=4)