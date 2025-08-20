import argparse
import pandas as pd
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import numpy as np



def regressor_metrics(y_test, y_pred):
    mae = mean_absolute_error(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    rmse = float(np.sqrt(mse))
    r2 = r2_score(y_test, y_pred)

    regression_metrics = {
        'r2': r2,
        'rmse': rmse,
        'mae': mae,
    }
    return regression_metrics


def calculate_metrics(args):
    # 从命令行参数中获取CSV文件路径
    df1_path = args.pred_path
    df2_path = args.test_path
    task = args.task
    data = args.data
    model_type = args.model
    output_csv = args.output_csv

    # 加载CSV文件
    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)
    
    # 获取预测值和真实值
    
    if task == 'reg':
        sort_list = ['logp',
                'pka_acidic',
                'cl-plasma',
                'logs',
                'logvdss',
                'pka_basic',
                't12',
                'logd',
                'mdck',
                'caco2',
                'ppb',
                'cl-int']


    df_reg = pd.DataFrame()
    for i in  sort_list:
        for m in range(3):
            task_reg_metrics = {}
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

            # 计算R^2（决定系数）
            reg_metrics = regressor_metrics(y_true, y_pred)
            task_reg_metrics[i] = reg_metrics
            df_task_reg = pd.DataFrame(task_reg_metrics).T
            df_reg = pd.concat([df_reg, df_task_reg])

    df_reg['task'] = df_reg.index
    grouped_reg_mean = df_reg.groupby('task').mean()
    grouped_reg_std = df_reg.groupby('task').std()
    grouped_reg_mean.reset_index(inplace=True)
    grouped_reg_std.reset_index(inplace=True)
    grouped_reg_mean = grouped_reg_mean.round(3)
    grouped_reg_std =grouped_reg_std.round(3)

    reg_columns = grouped_reg_mean.columns
    reg_mean_std = {}
    for col in reg_columns:
        reg_mean_std[col] = grouped_reg_mean[col].astype(str) + '±' + grouped_reg_std[col].astype(str)
    reg_mean_std['task'] = grouped_reg_mean['task']
    reg_mean_std = pd.DataFrame(reg_mean_std)
    reg_mean_std = reg_mean_std.set_index('task')
    reg_mean_std = reg_mean_std.loc[sort_list]
    reg_mean_std=reg_mean_std.reindex(columns=['r2','rmse','mae'])
    reg_mean_std.to_csv(output_csv)

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Calculate regression metrics")

    # 添加命令行参数
    parser.add_argument("--test_path", type=str, required=True, help="Path to the first CSV file")
    parser.add_argument("--pred_path", type=str, required=True, help="Path to the second CSV file")
    parser.add_argument("--task", type=str, required=True, help="task of data")
    parser.add_argument("--data", type=str, required=True, help="task of data")
    parser.add_argument("--model", type=str, required=True, help="model type")
    parser.add_argument("--output_csv", type=str, required=True, help="Path to save the output CSV file")
    # 解析命令行参数
    args = parser.parse_args()

    # 调用函数计算指标
    calculate_metrics(args)