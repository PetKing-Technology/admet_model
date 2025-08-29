# ADMET Property Prediction with Chemprop

This project provides a script to predict ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties of chemical compounds using a pre-trained Chemprop model.


## Installation & Environment Setup

Follow these steps to set up the necessary environment and dependencies.

1.  **Create and Activate Conda Environment**

    First, create a dedicated Conda environment for this project. We'll name it `admet_env` for clarity.

    ```bash
    # Create the environment with Python 3.11
    conda create -n admet_env python=3.11 -y

    # Activate the environment
    conda activate admet_env
    ```

2.  **Install Required Packages**

    Install the specific version of `chemprop` required to run the model.

    ```bash
    pip install chemprop==1.5.2
    pip install rdkit rdkit-contrib
    ```

## Usage

Follow these two steps to run a prediction.

### Step 1: Configure Model Path

You must first tell the script where to find the trained model file.

1.  Open the `predict.py` file in a text editor.
2.  Locate the `MODEL_PATHS` variable.
3.  Change the placeholder path to the **absolute path** of your model file.

    **Example:**

    ```python
    MODEL_PATHS = {
    'reg': "/path/to/your/model/reg_model_nordkit2d_2/fold_0",
    'cla': "/path/to/your/model/cla_model_nordkit2d/fold_0"
    }

    ```

### Step 2: Run Prediction

Once the path is correctly configured, execute the script from your terminal.

```bash
python predict.py
```

## Output

The script will process the input data and generate a file named `prediction_output.json` in the project's root directory. This file contains the ADMET predictions in JSON format.


### ⚡ 性能表现

#### 处理速度对比表

| 分子数量 | 总耗时 | 每秒处理 | 平均每分子 | 性能评估 |
|---------|--------|----------|------------|----------|
| 100     | 1.82s  | 54.9个   | 18.2ms     | ⭐⭐⭐⭐⭐ |
| 500     | 5.25s  | 95.2个   | 10.5ms     | ⭐⭐⭐⭐⭐ |
| 1,000   | 10.88s | 91.9个   | 10.9ms     | ⭐⭐⭐⭐⭐ |
| 5,000   | 118.85s| 42.1个   | 23.8ms     | ⭐⭐⭐⭐  |
| 10,000  | 413.03s| 24.2个   | 41.3ms     | ⭐⭐⭐   |


### 能计算的终点

```python
'Physicochemical_Property':['TPSA','nRot','nHA','nHD', 'MW', 'LogP', 'QED', 'SAscore', 'Lipinski', 'Pfizer', 'BMS', 'PAINS']
'cla': ['ames', 'lm-human', 'herg-10um', 'lm-mouse', 'f50', 'pgp_inh',
            'dili', 'herg', 'oatp1b1', 'cyp3a4-inh', 'cyp2d6-sub',
            'cyp3a4-sub', 'aggregators', 'oatp1b3', 'pgp_sub',
            'cyp2d6-inh', 'carcinogenicity', 'bbb'],
   'reg': ['logp', 'pka_acidic', 'cl-plasma', 'logs', 'logvdss',
            'pka_basic', 't12', 'logd', 'mdck', 'caco2', 'ppb', 'cl-int']
```

## Olaparib ADMET预测与实验结果对比分析

### 总结

本表格对比了Olaparib的ADMET预测结果与已发表的实验数据，评估预测模型的准确性。

| 参数 | 预测值 | 预测置信度 | 实验值/文献报告 | 预测准确性 | 备注 |
|------|---------|------------|----------------|------------|------|
| **药代动力学参数** |
| logP | 1.75 | Low-confidence | 2.35 (计算值) | ❌ 偏低 | 预测值比实际偏低约0.6个单位 |
| logD | 1.43 | Low-confidence | 1.49 (pH=7.4, 实验) | ✅ 接近 | 预测相对准确，误差约0.06 |
| logs (溶解性) | -3.33 | Low-confidence | 0.1 mg/mL = -3.6 log mg/mL | ✅ 接近 | 预测基本准确，实验确认为难溶性化合物 |
| pKa_acidic | 6.71 | Low-confidence | -1.16 (实验) | ❌ 完全错误 | 预测为弱酸性，实验显示强酸性 |
| pKa_basic | 6.78 | High-confidence | 12.07 (实验) | ❌ 严重偏低 | 预测严重低估碱性强度 |
| **药物分布参数** |
| PPB (血浆蛋白结合率) | 82.1% | High-confidence | 82-91% (实验) | ✅ 准确 | 预测值在实验范围内，高浓度时为82% |
| logVdss | -0.24 | High-confidence | V/F = 167L (实验) | ✅ 合理 | 预测的低分布容积与实验一致 |
| **代谢参数** |
| cl-plasma | 6.0 L/h | High-confidence | CL/F = 6.8 L/h (实验) | ✅ 接近 | 预测清除率与实验值接近 |
| cl-int | -1.40 | High-confidence | 主要通过CYP3A代谢 | ✅ 一致 | 与已知的CYP3A介导代谢一致 |
| t1/2 | 0.68 h | Low-confidence | 12 h (400mg剂量) | ❌ 严重偏低 | 预测半衰期过短，约为实际的1/18 |
| **渗透性参数** |
| Caco-2 | -5.36 | Low-confidence | BCS IV类化合物 | ✅ 一致 | 预测的低渗透性与BCS分类一致 |
| MDCK | -4.93 | Low-confidence | 低渗透性 | ✅ 一致 | 与实验观察的低渗透性一致 |
| **安全性参数** |
| Ames致突变性 | 0.27 (低风险) | Excellent | 无致突变性报告 | ✅ 准确 | 预测与安全性数据一致 |
| hERG (10μM) | 0.28 (低风险) | Excellent | 无心脏毒性报告 | ✅ 准确 | 临床未发现QT延长等心脏毒性 |
| hERG (常规) | 0.99 (高风险) | Poor | 无显著心脏毒性 | ❌ 假阳性 | 常规hERG预测过高估计风险 |
| DILI (肝毒性) | 0.72 (中等风险) | Poor | 罕见肝毒性报告 | ⚠️ 需谨慎 | 预测提示潜在风险，临床需监测 |
| 致癌性 | 0.08 (低风险) | Excellent | 无致癌性证据 | ✅ 准确 | 预测与长期安全性数据一致 |
| **转运蛋白相关** |
| P-gp底物 | 0.67 (中等可能) | Medium | P-gp底物 (实验确认) | ✅ 准确 | 与体外实验结果一致 |
| P-gp抑制剂 | 0.96 (高风险) | Poor | 非P-gp抑制剂 | ❌ 假阳性 | 预测过高估计抑制风险 |
| CYP3A4抑制剂 | 0.89 (高风险) | Poor | 弱CYP3A抑制剂 | ⚠️ 部分正确 | 预测方向正确但程度过高 |
| CYP2D6底物 | 0.24 (低可能) | Excellent | 主要非CYP2D6代谢 | ✅ 准确 | 与已知代谢途径一致 |
| **血脑屏障** |
| BBB透过性 | 0.35 (中等) | Medium | 有限脑部渗透 | ✅ 接近 | 与实验观察的有限脑暴露一致 |

## 预测模型表现评估

### 🎯 预测准确的参数 (8/25, 32%)
- 血浆蛋白结合率
- 分布容积 
- 清除率
- 溶解性
- 渗透性参数
- 基础安全性指标
- P-gp底物识别
- CYP2D6底物

### ⚠️ 需要改进的参数 (10/25, 40%)
- 脂水分配系数 (logP)
- pKa值预测
- 半衰期预测
- hERG风险评估
- 部分转运蛋白抑制预测

### ❌ 预测失败的参数 (7/25, 28%)
- pKa_acidic (完全错误)
- pKa_basic (严重偏差) 
- t1/2 (严重偏低)
- P-gp抑制剂 (假阳性)
- 常规hERG预测 (假阳性)

## 关键发现

1. **理化性质预测**: logD预测较准确，但pKa预测存在系统性错误，提示分子电离状态建模需要改进。

2. **药代动力学**: 清除率和分布容积预测较准确，但半衰期预测严重偏低，可能与模型未充分考虑非线性药代动力学有关。

3. **安全性评估**: 基础毒性预测(Ames、致癌性)表现良好，但hERG等特异性毒性预测存在假阳性问题。

4. **转运蛋白**: 底物识别能力较好，但抑制剂预测容易高估风险。

# 测试
# ADMET预测测评推荐药物列表

## 经典测试化合物
```
    compound_dict = {
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
```