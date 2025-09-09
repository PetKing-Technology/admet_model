import os
import json
import time
import pandas as pd
import chemprop
import math # 引入math库以处理可能的无效数字
from scripts.Process_runner import  ProcessRunner
from scripts.molecular_descriptors import MolecularDescriptors
from tqdm import tqdm  # 添加进度条支持
import logging  # 添加日志支持

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# --- 1. 配置区域 ---
MODEL_PATHS = {
    'reg': "/mnt/newdisk/fuli/ADMET/model_path/reg_model_nordkit2d_2/fold_0",
    'cla': "/mnt/newdisk/fuli/ADMET/model_path/cla_model_nordkit2d_3/fold_0"
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
   'cla': [
        "ames", "lm-human", "herg-10um", "lm-mouse",
        "f50", "pgp_inh", "dili", "oatp1b1", "cyp3a4-inh",
        "cyp2d6-sub", "H-HT", "Hepatotoxicity", "cyp3a4-sub",
        "aggregators", "oatp1b3", "pgp_sub", "cyp2d6-inh",
        "carcinogenicity", "bbb"
    ],
   'reg': ['logp', 'pka_acidic', 'cl-plasma', 'logs', 'logvdss',
            'pka_basic', 't12', 'logd', 'mdck', 'caco2', 'ppb', 'cl-int']
}


species_mapping = {
    'ames':'rat', 
    'lm-human':'human',
    'lm-mouse':'mouse',
    'herg-10um':'human', 
    'f50': 'human', 
    'pgp_inh': 'human',
    'pgp_sub':'human',
    'dili':'human', 
    'H-HT':'human',
    'Hepatotoxicity':'human',
    'oatp1b1':'human',
    'oatp1b3':'human',
    'cyp3a4-inh':'human', 
    'cyp2d6-sub':'human',
    'cyp3a4-sub':'human', 
    'cyp2d6-inh':'human', 
    'carcinogenicity':'human', 
    'cl-plasma':'human',
    'mdck':'dog', 
    'caco2':'human', 
    'ppb':'human', 
    'cl-int':'human',
    't12':'human', 
    'bbb':'human',
    'logvdss':'human',
    'aggregators':'common', 
    'logp':'common', 
    'pka_acidic':'common',  
    'logs':'common', 
    'pka_basic':'common', 'logd':'common', 
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

def validate_smiles_list(smiles_list):
    """
    验证SMILES列表的输入
    
    Args:
        smiles_list: 输入的SMILES列表
        
    Returns:
        tuple: (is_valid, error_message, cleaned_list)
    """
    if not isinstance(smiles_list, (list, tuple)):
        return False, "Input must be a list or tuple", []
    
    if len(smiles_list) == 0:
        return False, "SMILES list cannot be empty", []
    
    cleaned_list = []
    invalid_smiles = []
    
    for i, smiles in enumerate(smiles_list):
        if not isinstance(smiles, str):
            invalid_smiles.append(f"Index {i}: Not a string")
            continue
        
        if not smiles.strip():
            invalid_smiles.append(f"Index {i}: Empty string")
            continue
        
        # 基本的SMILES格式检查（简单验证）
        if len(smiles.strip()) < 3:
            invalid_smiles.append(f"Index {i}: Too short")
            continue
            
        cleaned_list.append(smiles.strip())
    
    if invalid_smiles:
        return False, f"Invalid SMILES found: {'; '.join(invalid_smiles[:5])}", []
    
    return True, "", cleaned_list

def chunk_list(lst, chunk_size):
    """
    将列表分割成指定大小的块
    
    Args:
        lst: 要分割的列表
        chunk_size: 每块的大小
        
    Returns:
        list: 分割后的块列表
    """
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

def classify_reg_confidence(uncertainty_value, threshold):
    """根据不确定度阈值判断回归任务的置信度"""

    output = []
    for i in uncertainty_value:
        if i > threshold:
            output.append("Low-confidence")
        else:
            output.append("High-confidence")
    return output


def classify_cla_confidence(probability):
    """根据预测概率判断分类任务的置信度"""
    if isinstance(probability, list):
        output = []
        for prob in probability:
            if not isinstance(prob, (float, int)):
                output.append("Invalid Input")
            elif prob < 0.3:
                output.append("excellent")
            elif prob < 0.7:
                output.append("medium")
            else:
                output.append("poor")
        return output
    else:
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
    支持单个值或值列表。
    """
    if isinstance(value, list):
        # 处理列表输入
        interpretations = []
        for val in value:
            interpretation = _get_single_interpretation(task_name, val)
            interpretations.append(interpretation)
        return interpretations
    else:
        # 处理单个值输入
        return _get_single_interpretation(task_name, value)

def _get_single_interpretation(task_name, value):
    """
    为单个值计算解释
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
            else: # value < 1
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

def calculate_molecular_descriptors(smiles_list):
    """
    计算分子描述符
    
    Args:
        smiles_list (list): SMILES字符串列表
        
    Returns:
        dict: 包含每个分子描述符的字典
    """

    try:
        calc = MolecularDescriptors()
        results = {}
        
        # 使用tqdm显示进度
        for smiles in tqdm(smiles_list, desc="Calculating molecular descriptors", unit="smiles"):
            try:
                result = calc.calculate_all_descriptors(smiles)
                if result:
                    results[smiles] = result
                else:
                    results[smiles] = None
            except Exception as e:
                logger.warning(f"Error calculating descriptors for {smiles}: {e}")
                results[smiles] = None
        
        return results
    except Exception as e:
        logger.error(f"Error calculating molecular descriptors: {e}")
        return {}

# --- 3. 核心预测函数 ---
# (此部分无变化)
def run_regression_predictions(smiles_list, smiles_file, pred_file):
    print("Running regression model predictions...")
    arguments = [
        '--test_path', smiles_file, '--preds_path', pred_file,
        '--checkpoint_dir', MODEL_PATHS['reg'], '--num_workers', '0',
        '--uncertainty_method', 'ensemble', '--gpu', '0'
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
        '--checkpoint_dir', MODEL_PATHS['cla'], '--num_workers', '0', '--gpu', '3'
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

    molecular_descriptors = calculate_molecular_descriptors(smiles_list)
    
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
            
            # 确保值是列表格式
            if hasattr(pred_val, '__iter__') and not isinstance(pred_val, str):
                pred_val = pred_val.tolist() if hasattr(pred_val, 'tolist') else list(pred_val)
            else:
                pred_val = [pred_val]
                
            if hasattr(un_val, '__iter__') and not isinstance(un_val, str):
                un_val = un_val.tolist() if hasattr(un_val, 'tolist') else list(un_val)
            else:
                un_val = [un_val]

            # 1. 创建基础的结果字典
            task_output = {
                "prediction": pred_val,
                "uncertainty": un_val,
                "confidence": classify_reg_confidence(un_val, UNCERTAINTY_THRESHOLDS.get(task))
            }
            task_output['species'] = species_mapping[task]
            # 2. 调用新函数获取解释，如果存在则添加到字典中
            interpretation = get_regression_interpretation(task, pred_val)
            if interpretation:
                task_output["interpretation"] = interpretation

            # 3. 将完整的任务结果存入最终输出
            final_json_output[smiles][task] = task_output

        # 处理分类结果 (无变化)
        for task in COLNAMES_DICT['cla']:
            prob_val = cla_results_df.loc[smiles, task]
            
            # 确保值是列表格式
            if hasattr(prob_val, '__iter__') and not isinstance(prob_val, str):
                prob_val = prob_val.tolist() if hasattr(prob_val, 'tolist') else list(prob_val)
            else:
                prob_val = [prob_val]
            
            final_json_output[smiles][task] = {
                "prediction": prob_val,
                "confidence": classify_cla_confidence(prob_val),
                'species': species_mapping[task]
            }
            
        
        # 添加分子描述符结果
        if molecular_descriptors and smiles in molecular_descriptors and molecular_descriptors[smiles]:
            final_json_output[smiles]['Physicochemical_Property'] = molecular_descriptors[smiles]
        else:
            final_json_output[smiles]['Physicochemical_Property'] = None
            
    return final_json_output

# --- 5. 脚本执行入口 ---
def load_smiles_from_file(file_path):
    """
    从文件加载SMILES列表
    
    Args:
        file_path (str): 文件路径，支持.txt和.csv格式
        
    Returns:
        list: SMILES字符串列表
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
            if 'smiles' in df.columns:
                return df['smiles'].dropna().tolist()
            else:
                return df.iloc[:, 0].dropna().tolist()  # 假设第一列是SMILES
        elif file_path.endswith('.txt'):
            with open(file_path, 'r') as f:
                return [line.strip() for line in f if line.strip()]
        else:
            raise ValueError("Unsupported file format. Please use .txt or .csv files.")
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {e}")

if __name__ == '__main__':
    # 确定SMILES列表
    compound_dict = {
    "Olaparib":"C1CC1C(=O)N1CCN(CC1)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C54)F",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "carbamazepine": "NC(=O)N1C2=CC=CC=C2C=C2C=CC=CC=C21",
    "digoxin": "CC1OC2CC3C4CCC5CC(O)CCC5(C)C4CCC3(C)C2(O)CC1OC1OC(C)C(O)C(O)C1O",
    "terfenadine": "CC(C)(C)C1=CC=C(C=C1)C(O)(CCCN2CCC(CC2)C(O)(C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5",
    "troglitazone": "CC1=C(C)C2=C(C=C1)OC(C)(COC3=CC=C(CC4SC(=O)NC4=O)C=C3)CC2",
    "cisapride": "COC1=CC(=CC(=C1OC)OC)C(=O)NC2CCN(CC2)CCCC3=CC=C(C=C3)F",
    "atorvastatin": "CC(C)C1=C(C(=C(N1CC(CC(=O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
    "imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    "omeprazole": "COC1=CC2=C(C=C1)N=C(N2)CS(=O)C3=NC4=C(N3)C=C(C=C4)OC",
    "clozapine": "CN1CCN(CC1)C2=NC3=CC=CC=C3NC4=C2C=C(C=C4)Cl",
    "paclitaxel": "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    "morphine": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
}

    test_smiles_list = list(compound_dict.values())
    # test_smiles_list = ["CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"]*100
    output = 'output_test.json'
    start_time = time.time()
    
    try:
        with suppress_stdout_stderr():
            final_predictions = main(test_smiles_list)
        end_time = time.time()
        
        
        # 保存结果
        with open(output, 'w') as f:
            json.dump(final_predictions, f, indent=4)
        print(f"mol: {len(test_smiles_list)} time: {end_time-start_time:.2f} seconds")

    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        print(f"Error: {e}")
        exit(1)