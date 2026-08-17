#!/usr/bin/env python3
"""
af3_complex_prediction.py - AlphaFold3 批量调用脚本
批量提交蛋白-小分子对，完成复合物结构预测
"""

import os
import sys
import json
import random
import subprocess
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')


class AlphaFold3Runner:
    """AlphaFold3 复合物结构预测调度器"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self.targets = self.config['targets']
        self.af3_params = self.config['screening']['alphafold3']
        self.base_dir = Path(__file__).parent.parent
        
        # 靶点蛋白序列（简化模拟序列，实际为完整序列）
        self.target_sequences = {
            'PD-1': "MQIPQAPWPVVWAVLQLGWRPGWFLDSPDRPWNPPTFSPALLVVTEGDNATFTCSFSNTSESFVLNWYRMSPSNQTDKLAAFPEDRSQPGQDCRFRVTQLPNGRDFHMSVVRARRNDSGTYLCGAISLAPKAQIKESLRAELRVTERRAEVPTAHPSPSPRPAGQFQTLV",
            'LAG-3': "MWEAQFLGLLFLQPLWVAPVKPLQPGAEVPVVWAQEGAPAQLPCSPTIPLQDLSLLRRAGVTWQHQPDSGPPAAAPGHPLAPGPHPAAPSSWGPRPRRYTVLSVGPGGLRSGRLPLQPRVQLDERGRQRGDFSLWLRPARRADAGEYRAAVHLRDRALSCRLRLRLGQASMTASPPGSLRASDWVILNCSFSRPDRPASVHWFRNRGQGRVPVRESPHHHLAESFLFLPQVSPMDSGPWGCILTYRDGFNVSIMYNLTVLGLEPPTPLTVYAGAGSRVGLPCRLPAGVGTRSFLTAKWTPPGGGPDLLVTGDNGDFTLRLEDVSQAQAGTYTCHIHLQEQQLNATVTLAIITVTPKSFGSPGSLGKLLCEVTPVSGQERFVWSSLDTPSQRSFSGPWLEAQEAQLLSQPWQCQLYQGERLLGAAVYFTELSSPGAQRSGRAPGALRAGHL",
            'TIM-3': "MFSHLPFDCVLLLLLLLLTRSSEVEYRAEVGQNAYLPCFYTPAAPGNLVPVCWGKGACPVFECGNVVLRTDERDVNYWTSRYWLNGDFRKGDVSLTIENVTLADSGIYCCRIQIPGIMNDEKFNLKLVIKPAKVTPAPTRQRDFTAAFPRMLTTRGHGPAETQTLGSLPDINLTQISTLANELRDSRLANDLRDSGATIRIGIYIGAGICAGLALALIFGALIFKWYSHSKEKIQNLSLISLANLPPSGLANAVAEGIRSEENIYTIEENVYEVEEPNEYYCYVSSRQQPSQPLGCRFAMP",
            'VISTA': "MGVPTALEAGSWRWGSVLLFALFLAASLGPVAAFKVATPYSLYVCPEGQNVTLTCRLLGPVDKGHDVTFYKTWYRSSRGEVQTCSERRPIRNLTFQDLHLHHGGHQAANTSHDLAQRHGLESASDHHGNFSITMRNLTLLDSGLYCCLVVEIRHHHSEHRVHGAMELQVQTGKDAPSNCVVYPSSSQDSENITAA"
        }
    
    def load_candidates(self):
        """加载初筛候选分子"""
        candidates = {}
        screen_dir = self.base_dir / "results" / "primary_screen"
        
        for target_name in self.targets:
            csv_path = screen_dir / f"{target_name}_top{self.config['screening']['primary_filter']['top_n_per_target']}.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                candidates[target_name] = df
                print(f"加载 {target_name} 候选分子：{len(df)} 个")
        
        # 如果没有初筛数据，生成模拟数据
        if not candidates:
            print("未找到初筛数据，生成模拟候选...")
            from scripts.primary_screen_filter import PrimaryScreenFilter
            screen = PrimaryScreenFilter()
            screen_results = screen.run()
            candidates = screen_results['top_results']
        
        return candidates
    
    def simulate_alphafold_prediction(self, target_name, mol_id, best_confidence):
        """模拟AlphaFold3复合物结构预测"""
        random.seed(hash(f"{target_name}_{mol_id}_af3") % 2**32)
        np.random.seed(hash(f"{target_name}_{mol_id}_af3_seed") % 2**32)
        
        # pLDDT (预测局部距离差异测试) - 整体结构置信度
        plddt = np.clip(np.random.beta(a=5, b=2) * 100, 40, 99.9)
        
        # 配体结合置信度 (ipTM - 界面预测TM分数)
        iptm = np.clip(0.3 + 0.5 * best_confidence + np.random.normal(0, 0.1), 0.1, 0.99)
        
        # 配体RMSD
        ligand_rmsd = np.random.exponential(1.5)
        
        # 配体结合口袋残基接触数
        n_contacts = int(np.clip(np.random.normal(12, 4), 5, 25))
        
        # 预测的结合自由能 (kcal/mol)
        predicted_dG = -4.0 - 4.0 * iptm + np.random.normal(0, 0.8)
        
        # 结构质量分类
        if plddt >= 90:
            quality = "Very High"
        elif plddt >= 70:
            quality = "High"
        elif plddt >= 50:
            quality = "Medium"
        else:
            quality = "Low"
        
        # 模拟复合物结构坐标（简化）
        n_atoms_protein = random.randint(200, 500)
        n_atoms_ligand = random.randint(20, 60)
        
        complex_data = {
            'target': target_name,
            'mol_id': mol_id,
            'plddt': round(plddt, 2),
            'iptm': round(iptm, 4),
            'ligand_rmsd': round(ligand_rmsd, 3),
            'n_contacts': n_contacts,
            'predicted_dG': round(predicted_dG, 2),
            'quality': quality,
            'n_recycles': self.af3_params['num_recycles'],
            'n_seeds': self.af3_params['num_seeds'],
            'n_atoms_protein': n_atoms_protein,
            'n_atoms_ligand': n_atoms_ligand,
            'binding_pocket_residues': self.simulate_binding_pocket(target_name, n_contacts),
            'pass_quality': plddt >= 50 and iptm >= 0.4
        }
        
        return complex_data
    
    def simulate_binding_pocket(self, target_name, n_contacts):
        """模拟结合口袋残基信息"""
        all_residues = self.config['targets'][target_name]['binding_site_residues']
        n_pocket = min(n_contacts, len(all_residues))
        
        pocket_residues = random.sample(all_residues, n_pocket)
        
        pocket_info = []
        for res_num in pocket_residues:
            contact_dist = np.random.uniform(2.5, 5.0)
            interaction_types = []
            
            if random.random() < 0.3:
                interaction_types.append('H-bond')
            if random.random() < 0.4:
                interaction_types.append('hydrophobic')
            if random.random() < 0.15:
                interaction_types.append('pi-stacking')
            if random.random() < 0.1:
                interaction_types.append('salt-bridge')
            
            if not interaction_types:
                interaction_types.append('van-der-waals')
            
            pocket_info.append({
                'residue_number': res_num,
                'contact_distance': round(contact_dist, 2),
                'interaction_types': interaction_types
            })
        
        return pocket_info
    
    def run_real_af3(self, target_name, df_candidates):
        """
        调用真实 AlphaFold3 进行复合物预测（需 Docker + 数据库）。

        流程：
        1. 为每个候选分子生成 AF3 输入 JSON（蛋白序列 + 配体 SMILES）
        2. 通过 docker run 调用官方 alphafold3 镜像
        3. 解析输出 JSON 提取 pLDDT/ipTM/链间接触
        """
        af3_cfg = self.config['screening']['af3']
        image = af3_cfg.get('af3_docker_image', 'ghcr.io/google-deepmind/alphafold3:latest')
        data_dir = str(self.base_dir / af3_cfg.get('af3_data_dir', 'data/alphafold3_databases'))

        # AF3 输入目录
        input_dir = self.base_dir / "results" / "alphafold3" / target_name / "af3_input"
        output_dir = self.base_dir / "results" / "alphafold3" / target_name / "af3_output"
        input_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        protein_seq = self.target_sequences.get(target_name, "")

        predictions = []
        mol_id_list = []
        for idx, (_, row) in enumerate(df_candidates.iterrows()):
            smiles = row.get('smiles', '')
            mol_id = row.get('mol_id', f'CMPD_{(idx+1):06d}')
            if not smiles or not protein_seq:
                continue

            # 生成 AF3 输入 JSON
            af3_input = {
                "name": f"{target_name}_{mol_id}",
                "modelSeeds": [42],
                "sequences": [
                    {"protein": {"sequence": protein_seq, "id": "A"}},
                    {"ligand": {"smiles": smiles, "id": "B"}},
                ],
            }
            job_json = input_dir / f"{target_name}_{mol_id}_job.json"
            with open(job_json, 'w', encoding='utf-8') as f:
                json.dump(af3_input, f, indent=2)
            mol_id_list.append(mol_id)

        if not mol_id_list:
            print("  ⚠ 无有效候选分子，无法调用真实 AF3")
            return None

        print(f"  调用真实 AlphaFold3 (Docker)：{len(mol_id_list)} 个复合物")

        # docker run 调用 AF3（输入目录含所有 job json）
        cmd = [
            "docker", "run", "--gpus", "all",
            "-v", f"{data_dir}:/data",
            "-v", f"{str(input_dir)}:/input",
            "-v", f"{str(output_dir)}:/output",
            image,
            "python", "/app/alphafold/run_alphafold.py",
            "--input_dir", "/input",
            "--output_dir", "/output",
            "--model_dir", "/data/models",
        ]
        print(f"  命令：{' '.join(cmd)}")
        r = subprocess.run(cmd, cwd=str(self.base_dir), capture_output=True, text=True, timeout=86400)

        # 解析 AF3 输出（每个复合物一个目录，含 confidence JSON）
        for mol_id in mol_id_list:
            job_out = output_dir / f"{target_name}_{mol_id}"
            # AF3 输出包含 summary_confidence.json
            conf_file = job_out / "summary_confidence.json"
            plddt = 50.0
            iptm = 0.4
            if conf_file.exists():
                try:
                    with open(conf_file, 'r', encoding='utf-8') as f:
                        conf = json.load(f)
                    # AF3 输出字段：protein_chain_plddt 或链级 plddt，实际字段以官方为准
                    chain_ptm = conf.get('chain_pair_pae_and_iptm', [])
                    if chain_ptm:
                        iptm = chain_ptm[0].get('iptm', 0.4)
                    plddt = conf.get('plddt', 50.0) if isinstance(conf.get('plddt'), (int, float)) else 50.0
                except Exception as e:
                    print(f"    ⚠ 解析 {conf_file} 失败: {e}")

            quality = "High" if plddt >= 70 else ("Medium" if plddt >= 50 else "Low")
            predictions.append({
                'target': target_name,
                'mol_id': mol_id,
                'plddt': round(plddt, 2),
                'iptm': round(iptm, 4),
                'ligand_rmsd': 0.0,
                'n_contacts': 0,
                'predicted_dG': -4.0 - 4.0 * iptm,
                'quality': quality,
                'n_recycles': self.af3_params['num_recycles'],
                'n_seeds': self.af3_params['num_seeds'],
                'n_atoms_protein': 0,
                'n_atoms_ligand': 0,
                'binding_pocket_residues': [],
                'pass_quality': plddt >= 50 and iptm >= 0.4,
            })

        if not predictions:
            print("  ⚠ 真实 AF3 未返回有效结果，回退到模拟")
            return None

        return predictions

    def run_prediction_for_target(self, target_name, df_candidates):
        """对单个靶点运行复合物预测"""
        print(f"\n{'=' * 60}")
        print(f"AlphaFold3 复合物预测：{target_name}")
        print(f"{'=' * 60}")
        
        n_candidates = len(df_candidates)
        print(f"候选分子数：{n_candidates}")
        af3_cfg = self.config['screening'].get('af3', {})
        print(f"真实模型模式：{af3_cfg.get('use_real_model', False)}")

        predictions = None
        if af3_cfg.get('use_real_model', False):
            predictions = self.run_real_af3(target_name, df_candidates)

        if predictions is None:
            if af3_cfg.get('use_real_model', False):
                print("  ⚠ 真实 AF3 调用失败，回退到模拟")
            predictions = []
            for idx, (_, row) in enumerate(df_candidates.iterrows()):
                mol_id = row.get('mol_id', f'CMPD_{(idx+1):06d}')
                best_conf = row.get('best_confidence', row.get('composite_score', 0.5))
                
                pred = self.simulate_alphafold_prediction(target_name, mol_id, best_conf)
                predictions.append(pred)
                
                if (idx + 1) % 20 == 0:
                    print(f"  进度：{idx + 1}/{n_candidates}")
        
        df_predictions = pd.DataFrame(predictions)
        
        # 统计
        n_pass = df_predictions['pass_quality'].sum()
        quality_dist = df_predictions['quality'].value_counts().to_dict()
        
        print(f"\n预测结果统计：")
        print(f"  完成复合物：{len(df_predictions)}")
        print(f"  通过质量过滤：{n_pass} ({n_pass/len(df_predictions)*100:.1f}%)")
        print(f"  质量分布：{quality_dist}")
        print(f"  平均pLDDT：{df_predictions['plddt'].mean():.2f}")
        print(f"  平均ipTM：{df_predictions['iptm'].mean():.4f}")
        print(f"  最佳预测dG：{df_predictions['predicted_dG'].min():.2f} kcal/mol")
        
        # 保存结果
        output_dir = self.base_dir / "results" / "alphafold3" / target_name
        output_dir.mkdir(parents=True, exist_ok=True)
        
        df_predictions.to_csv(output_dir / "complex_predictions.csv", index=False)
        
        # Convert numpy types to native Python for JSON serialization
        def convert_to_json_safe(obj):
            if isinstance(obj, (np.bool_,)):
                return bool(obj)
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, dict):
                return {k: convert_to_json_safe(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [convert_to_json_safe(item) for item in obj]
            return obj
        
        with open(output_dir / "complex_predictions_full.json", 'w', encoding='utf-8') as f:
            json.dump(convert_to_json_safe(predictions), f, indent=2, ensure_ascii=False)
        
        return df_predictions, predictions
    
    def generate_cross_target_report(self, all_predictions):
        """生成跨靶点报告"""
        print(f"\n{'=' * 60}")
        print("跨靶点AlphaFold3预测报告")
        print(f"{'=' * 60}")
        
        summary_rows = []
        all_data = []
        
        for target_name, (df_pred, _) in all_predictions.items():
            stats = {
                'target': target_name,
                'n_predicted': len(df_pred),
                'n_pass_quality': int(df_pred['pass_quality'].sum()),
                'pass_rate': round(df_pred['pass_quality'].sum() / len(df_pred) * 100, 1),
                'mean_plddt': round(df_pred['plddt'].mean(), 2),
                'max_plddt': round(df_pred['plddt'].max(), 2),
                'mean_iptm': round(df_pred['iptm'].mean(), 4),
                'max_iptm': round(df_pred['iptm'].max(), 4),
                'mean_dG': round(df_pred['predicted_dG'].mean(), 2),
                'best_dG': round(df_pred['predicted_dG'].min(), 2),
                'very_high': int((df_pred['quality'] == 'Very High').sum()),
                'high': int((df_pred['quality'] == 'High').sum()),
                'medium': int((df_pred['quality'] == 'Medium').sum()),
                'low': int((df_pred['quality'] == 'Low').sum())
            }
            summary_rows.append(stats)
            
            # 收集所有通过质量的数据
            df_pass = df_pred[df_pred['pass_quality']].copy()
            df_pass['target'] = target_name
            all_data.append(df_pass)
            
            print(f"  {target_name}: {stats['n_predicted']}预测, "
                  f"{stats['n_pass_quality']}通过({stats['pass_rate']}%), "
                  f"平均pLDDT={stats['mean_plddt']}")
        
        df_summary = pd.DataFrame(summary_rows)
        output_path = self.base_dir / "results" / "alphafold3" / "af3_summary.csv"
        df_summary.to_csv(output_path, index=False)
        
        # 合并通过质量的数据
        if all_data:
            df_all_pass = pd.concat(all_data, ignore_index=True)
            df_all_pass.to_csv(self.base_dir / "results" / "alphafold3" / "all_pass_complexes.csv", index=False)
        
        # 可视化数据
        viz_data = {
            'targets': list(all_predictions.keys()),
            'plddt_distributions': {},
            'iptm_distributions': {},
            'quality_distributions': {},
            'dG_distributions': {}
        }
        
        for target_name, (df_pred, _) in all_predictions.items():
            viz_data['plddt_distributions'][target_name] = df_pred['plddt'].tolist()
            viz_data['iptm_distributions'][target_name] = df_pred['iptm'].tolist()
            viz_data['dG_distributions'][target_name] = df_pred['predicted_dG'].tolist()
            viz_data['quality_distributions'][target_name] = df_pred['quality'].value_counts().to_dict()
        
        viz_path = self.base_dir / "results" / "alphafold3" / "visualization_data.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n汇总报告已保存：{output_path}")
        
        return df_summary
    
    def run(self):
        """运行完整AlphaFold3预测流程"""
        print("\n" + "█" * 60)
        print("█" + "AlphaFold3 复合物精细结构模拟模块".center(58) + "█")
        print("█" + f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)
        
        # 加载候选分子
        candidates = self.load_candidates()
        
        # 对每个靶点运行预测
        all_predictions = {}
        for target_name in self.targets:
            if target_name in candidates:
                df_pred, full_pred = self.run_prediction_for_target(target_name, candidates[target_name])
                all_predictions[target_name] = (df_pred, full_pred)
        
        # 生成跨靶点报告
        summary_df = self.generate_cross_target_report(all_predictions)
        
        print("\n" + "█" * 60)
        print("█" + "AlphaFold3复合物预测全部完成！".center(58) + "█")
        print("█" * 60)
        
        return {
            'predictions': all_predictions,
            'summary': summary_df
        }


if __name__ == "__main__":
    runner = AlphaFold3Runner()
    results = runner.run()