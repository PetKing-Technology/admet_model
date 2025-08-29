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

