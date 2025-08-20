import argparse
import pandas as pd
from sklearn.metrics import accuracy_score 
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_score
from sklearn.metrics import r2_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import recall_score
from sklearn.metrics import matthews_corrcoef
import numpy as np


def classifier_metrics(y_true, y_pred_s):
    # 模型的苹评估
    y_pred = []
    y_test = []
    y_pred_p = []
    for s, t in zip(y_pred_s, y_true):
        try:
            s = float(s)
            if s >= 0.5:
                y_pred.append(1)
            else:
                y_pred.append(0)
            y_test.append(t)
            y_pred_p.append(s)
        except ValueError:
            print("无法将字符串转换为浮点数：无效的输入")
    acc = accuracy_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred_p)
    TN, FP, FN, TP = confusion_matrix(y_test, y_pred).ravel()
    SE = TP / (TP + FN) # 敏感度是实际为正类的样本中，被模型正确预测为正类的比例
    SP = TN / (TN + FP) # 特异性是实际为负类的样本中，被模型正确预测为负类的比例。
    MCC = matthews_corrcoef(y_test, y_pred)
    f1_label1 = f1_score(y_test, y_pred, pos_label=1)
    precision_label1 = precision_score(y_test, y_pred, pos_label=1)
    f1_label0 = f1_score(y_test, y_pred, pos_label=0)
    precision_label0 = precision_score(y_test, y_pred, pos_label=0)
    classification_metrics = {
       'AUC':round(roc_auc,4),
       'ACC': round(acc, 4),
       'SP':round(SP,4),
       'SE':round(SE, 4),
       'MCC':round(MCC,4),
       'f1_label1':round(f1_label1,4),
       'f1_label0':round(f1_label0,4),
       'precision_label1':round(precision_label1,4),
       'precision_label0':round(precision_label0,4)
    }
    return classification_metrics


def calculate_metrics(args):
    # 从命令行参数中获取CSV文件路径
    df1_path = args.pred_path
    df2_path = args.test_path
    task = args.task
    data = args.data
    model_type = args.model
    output_csv = args.output_csv

    if task == 'cla':
        sort_list = ['ames',
                'lm-human',
                'herg-10um',
                'lm-mouse',
                'f50',
                'pgp_inh',
                'dili',
                'herg',
                'oatp1b1',
                'cyp3a4-inh',
                'cyp2d6-sub',
                'cyp3a4-sub',
                'aggregators',
                'oatp1b3',
                'pgp_sub',
                'cyp2d6-inh',
                'carcinogenicity',
                'bbb']


    # 加载CSV文件
    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)
    
    task_cla_metrics = {}
    # 获取预测值和真实值
    df_cla = pd.DataFrame()
    for i in sort_list:
        for m in range(3):
            if m == 0:
                col_name = i
            if m == 1:
                col_name = f'{i}_model_0'
            if m == 2:
                col_name = f'{i}_model_1'
            index = df2[i].dropna().index
            assert list(df1.loc[index,"smiles"]) == list(df2.loc[index,"smiles"])
            
            y_pred = df1.loc[index, col_name].values
            y_true = df2.loc[index, i].values
            task_cla_metrics[i] = classifier_metrics(y_true, y_pred)
            df_task_cla = pd.DataFrame(task_cla_metrics).T
            df_cla = pd.concat([df_cla, df_task_cla])

    df_cla['task'] = df_cla.index
    grouped_cla_mean = df_cla.groupby('task').mean()
    grouped_cla_std = df_cla.groupby('task').std()
    grouped_cla_mean.reset_index(inplace=True)
    grouped_cla_std.reset_index(inplace=True)
    grouped_cla_mean =grouped_cla_mean.round(3)
    grouped_cla_std =grouped_cla_std.round(3)

    cla_columns = grouped_cla_mean.columns
    cla_mean_std = {}
    for col in cla_columns:
        cla_mean_std[col] = grouped_cla_mean[col].astype(str) + '±' + grouped_cla_std[col].astype(str)
    cla_mean_std['task'] = grouped_cla_mean['task']
    cla_mean_std = pd.DataFrame(cla_mean_std)
    cla_mean_std = cla_mean_std.set_index('task')
    cla_mean_std = cla_mean_std.loc[sort_list]
    cla_mean_std=cla_mean_std.reindex(columns=[
        "AUC", "ACC","SP","SE","MCC","f1_label1","f1_label0", "precision_label1", "precision_label0"
    ])
    cla_mean_std.to_csv(output_csv)



if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Calculate regression metrics")

    # 添加命令行参数
    parser.add_argument("--test_path", type=str, required=True, help="Path to the first CSV file")
    parser.add_argument("--pred_path", type=str, required=True, help="Path to the second CSV file")
    parser.add_argument("--task", type=str, required=True, help="Path to the second CSV file")
    parser.add_argument("--data", type=str, required=True, help="Path to the second CSV file")
    parser.add_argument("--model", type=str, required=True, help="model type")
    parser.add_argument("--output_csv", type=str, required=True, help="Path to save the output CSV file")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用函数计算指标
    calculate_metrics(args)