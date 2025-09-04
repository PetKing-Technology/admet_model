"""
分子描述符和药物性质计算模块

该模块提供了计算以下分子描述符和药物性质的功能：
- TPSA（极性表面积）
- nRot（可旋转键数）
- nHA（氢键受体数）
- nHD（氢键供体数）
- Lipinski Rule（口服药物规则）
- Pfizer Rule（毒性预测规则）
- BMS Rule
- PAINS（假阳性报警结构）- 如果可用的话
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import os
import sys

# 检查可选模块的可用性
QED_AVAILABLE = False
SA_SCORE_AVAILABLE = False
PAINS_AVAILABLE = False

try:
    from rdkit.Chem import QED
    QED_AVAILABLE = True
except ImportError:
    print("Warning: QED module not available. Install rdkit-contrib for QED support.")


# 尝试导入SA Score和PAINS
try:
    # 添加RDKit Contrib路径
    import site
    conda_prefix = os.path.dirname(site.getsitepackages()[0])
    if conda_prefix:
        contrib_path = os.path.join(conda_prefix, 'site-packages', 'rdkit', 'Contrib')
        if os.path.exists(contrib_path):
            sys.path.append(contrib_path)
            print(f"Added RDKit Contrib path: {contrib_path}")
        else:
            print(f"RDKit Contrib path not found: {contrib_path}")
    
    # 尝试导入SA Score
    try:
        from SA_Score import sascorer
        SA_SCORE_AVAILABLE = True
        print("SA Score module imported successfully")
    except ImportError as e:
        print(f"Warning: SA Score module not available: {e}")
    
    # 不再需要NP Score，使用标准RDKit PAINS过滤器
    PAINS_AVAILABLE = True
    print("PAINS filters available via standard RDKit")
        
except Exception as e:
    print(f"Warning: Error setting up RDKit Contrib modules: {e}")

class MolecularDescriptors:
    """分子描述符计算类"""
    
    def __init__(self):
        """初始化分子描述符计算器"""
        self.pains_catalog = None
        if PAINS_AVAILABLE:
            self._init_pains_catalog()
    
    def _init_pains_catalog(self):
        """初始化PAINS过滤器目录"""
        try:
            params_pains = FilterCatalogParams()
            # 使用完整的PAINS目录
            params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
            self.pains_catalog = FilterCatalog(params_pains)
            print("PAINS catalog initialized successfully")
        except Exception as e:
            print(f"Warning: Could not initialize PAINS catalog: {e}")
            self.pains_catalog = None
    
    def calculate_all_descriptors(self, smiles):
        """
        计算所有分子描述符
        
        Args:
            smiles (str): SMILES字符串
            
        Returns:
            dict: 包含所有描述符的字典
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        results = {}
        
        # 基础分子描述符
        results['TPSA'] = self.calculate_tpsa(mol)
        results['nRot'] = self.calculate_nrot(mol)
        results['nHA'] = self.calculate_nha(mol)
        results['nHD'] = self.calculate_nhd(mol)
        results['MW'] = self.calculate_mw(mol)

        
        # 药物相似性评分（如果可用）
        if QED_AVAILABLE:
            results['QED'] = self.calculate_qed(mol)
        else:
            results['QED'] = None
            
        if SA_SCORE_AVAILABLE:
            results['SAscore'] = self.calculate_sascore(mol)
        else:
            results['SAscore'] = None
        
        # 药物规则
        results['Lipinski'] = self.check_lipinski_rule(mol)
        results['Pfizer'] = self.check_pfizer_rule(mol)
        results['BMS'] = self.check_bms_rule(mol)
        
        # PAINS检查
        results['PAINS'] = self.check_pains(mol)['PAINS']
        
        return results
    
    def calculate_tpsa(self, mol):
        """计算拓扑极性表面积 (TPSA)"""
        try:
            return Descriptors.TPSA(mol)
        except:
            return None
    
    def calculate_nrot(self, mol):
        """计算可旋转键数"""
        try:
            return Descriptors.NumRotatableBonds(mol)
        except:
            return None
    
    def calculate_nha(self, mol):
        """计算氢键受体数"""
        try:
            return Descriptors.NumHAcceptors(mol)
        except:
            return None
    
    def calculate_nhd(self, mol):
        """计算氢键供体数"""
        try:
            return Descriptors.NumHDonors(mol)
        except:
            return None
    
    def calculate_mw(self, mol):
        """计算分子量"""
        try:
            return Descriptors.MolWt(mol)
        except:
            return None
    
    
    def calculate_qed(self, mol):
        """计算药物相似性评分 (QED)"""
        if not QED_AVAILABLE:
            return None
        try:
            return QED.default(mol)
        except:
            return None
    
    def calculate_sascore(self, mol):
        """计算合成可行性评分 (SA Score)"""
        if not SA_SCORE_AVAILABLE:
            return None
        try:
            return sascorer.calculateScore(mol)
        except:
            return None
    
    def check_lipinski_rule(self, mol):
        """
        检查Lipinski规则（口服药物规则）
        
        Returns:
            dict: 包含各项规则检查结果的字典
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            nha = Descriptors.NumHAcceptors(mol)
            nhd = Descriptors.NumHDonors(mol)
            
            rules = {
                'MW': {'value': mw, 'pass': mw <= 500, 'rule': 'MW ≤ 500'},
                'LogP': {'value': logp, 'pass': logp <= 5, 'rule': 'LogP ≤ 5'},
                'HBA': {'value': nha, 'pass': nha <= 10, 'rule': 'HBA ≤ 10'},
                'HBD': {'value': nhd, 'pass': nhd <= 5, 'rule': 'HBD ≤ 5'}
            }
            
            # 计算通过的规则数
            passed_rules = sum(1 for rule in rules.values() if rule['pass'])
            rules['total_passed'] = passed_rules
            rules['overall_pass'] = passed_rules >= 3  # 至少通过3条规则
            
            return rules
        except:
            return None
    
    def check_pfizer_rule(self, mol):
        """
        检查Pfizer规则（毒性预测规则）
        
        Returns:
            dict: 包含各项规则检查结果的字典
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            nha = Descriptors.NumHAcceptors(mol)
            nhd = Descriptors.NumHDonors(mol)
            nrot = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            
            rules = {
                'MW': {'value': mw, 'pass': mw <= 400, 'rule': 'MW ≤ 400'},
                'LogP': {'value': logp, 'pass': logp <= 4.15, 'rule': 'LogP ≤ 4.15'},
                'HBA': {'value': nha, 'pass': nha <= 8, 'rule': 'HBA ≤ 8'},
                'HBD': {'value': nhd, 'pass': nhd <= 5, 'rule': 'HBD ≤ 5'},
                'nRot': {'value': nrot, 'pass': nrot <= 10, 'rule': 'nRot ≤ 10'},
                'TPSA': {'value': tpsa, 'pass': tpsa <= 90, 'rule': 'TPSA ≤ 90'}
            }
            
            # 计算通过的规则数
            passed_rules = sum(1 for rule in rules.values() if rule['pass'])
            rules['total_passed'] = passed_rules
            rules['overall_pass'] = passed_rules >= 4  # 至少通过4条规则
            
            return rules
        except:
            return None
    
    def check_bms_rule(self, mol):
        """
        检查BMS规则
        
        Returns:
            dict: 包含各项规则检查结果的字典
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            nha = Descriptors.NumHAcceptors(mol)
            nhd = Descriptors.NumHDonors(mol)
            nrot = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            
            rules = {
                'MW': {'value': mw, 'pass': mw <= 400, 'rule': 'MW ≤ 400'},
                'LogP': {'value': logp, 'pass': logp <= 4.0, 'rule': 'LogP ≤ 4.0'},
                'HBA': {'value': nha, 'pass': nha <= 8, 'rule': 'HBA ≤ 8'},
                'HBD': {'value': nhd, 'pass': nhd <= 3, 'rule': 'HBD ≤ 3'},
                'nRot': {'value': nrot, 'pass': nrot <= 8, 'rule': 'nRot ≤ 8'},
                'TPSA': {'value': tpsa, 'pass': tpsa <= 90, 'rule': 'TPSA ≤ 90'}
            }
            
            # 计算通过的规则数
            passed_rules = sum(1 for rule in rules.values() if rule['pass'])
            rules['total_passed'] = passed_rules
            rules['overall_pass'] = passed_rules >= 4  # 至少通过4条规则
            
            return rules
        except:
            return None
    
    def check_pains(self, mol):
        """
        检查PAINS（假阳性报警结构）
        
        Returns:
            dict: 包含PAINS检查结果的字典
        """
        if not PAINS_AVAILABLE:
            return {'available': False, 'message': 'PAINS module not available'}
        
        try:
            # 使用RDKit PAINS过滤器进行检查
            matches = []
            if self.pains_catalog:
                flag = self.pains_catalog.HasMatch(mol)

            
            return {
                'PAINS': flag,
            }
        except Exception as e:
            return {'PAINS': f'Error during PAINS check: {e}'}

def calculate_descriptors_batch(smiles_list):
    """
    批量计算分子描述符
    
    Args:
        smiles_list (list): SMILES字符串列表
        
    Returns:
        list: 包含每个分子描述符的列表
    """
    calculator = MolecularDescriptors()
    results = []
    
    for smiles in smiles_list:
        result = calculator.calculate_all_descriptors(smiles)
        results.append({
            'smiles': smiles,
            'descriptors': result
        })
    
    return results

# 使用示例
if __name__ == "__main__":
    # 测试SMILES
    test_smiles = "CC(C)OC(=O)CC(=O)CSc1nc2c(cc1C#N)CCC2"
    
    # 创建计算器实例
    calc = MolecularDescriptors()
    
    # 计算单个分子的描述符
    result = calc.calculate_all_descriptors(test_smiles)
    
    if result:
        print("分子描述符计算结果:")
        print(f"SMILES: {test_smiles}")
        print(f"分子量: {result['MW']:.2f}")
        print(f"LogP: {result['LogP']:.2f}")
        print(f"TPSA: {result['TPSA']:.2f}")
        print(f"可旋转键数: {result['nRot']}")
        print(f"氢键受体数: {result['nHA']}")
        print(f"氢键供体数: {result['nHD']}")
        
        if result['QED'] is not None:
            print(f"QED: {result['QED']:.3f}")
        else:
            print("QED: 不可用")
            
        if result['SAscore'] is not None:
            print(f"SA Score: {result['SAscore']:.2f}")
        else:
            print("SA Score: 不可用")
        
        print("\nLipinski规则检查:")
        for rule, data in result['Lipinski'].items():
            if rule not in ['total_passed', 'overall_pass']:
                status = "✓" if data['pass'] else "✗"
                print(f"  {rule}: {data['value']:.2f} {status} ({data['rule']})")
        print(f"  总体: {'通过' if result['Lipinski']['overall_pass'] else '不通过'}")
        
        print("\nPAINS检查:")
        pains_result = result['PAINS']
        if pains_result['PAINS'] == True:
            print("  发现PAINS结构")
        else:
            print("  未发现PAINS结构")
    else:
        print("无法解析SMILES字符串")
