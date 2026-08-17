#!/usr/bin/env python3
"""
tme_simulator.py - 肿瘤微环境 (TME) 模拟模块
模拟酸性缺氧环境对蛋白-配体结合的影响，输出TME修正因子
"""

import os
import sys
import json
import random
import numpy as np
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')


# ========================================
# 肿瘤微环境 (TME) 参数
# ========================================
"""
肿瘤微环境 vs 正常组织环境对比：

┌─────────────────────┬──────────────┬─────────────┬──────────────────────────┐
│ 环境参数             │ 正常组织      │ 肿瘤微环境   │ 对蛋白-配体结合的影响      │
├─────────────────────┼──────────────┼─────────────┼──────────────────────────┤
│ pH                   │ 7.35 - 7.45  │ 6.2 - 6.8   │ 质子化状态改变             │
│ 氧浓度 (pO2)         │ 40-60 mmHg   │ 0-20 mmHg   │ 缺氧诱导因子(HIF)激活      │
│ 温度                 │ 37°C         │ 35-38°C     │ 动力学速率变化             │
│ 乳酸浓度             │ <2 mM        │ 10-40 mM    │ 竞争结合位点               │
│ 活性氧 (ROS)         │ 低           │ 高           │ 氧化修饰                   │
│ 免疫抑制因子         │ 低           │ 高           │ 竞争/协同效应              │
│ (IL-10, TGF-β, VEGF) │              │             │                           │
│ ATP 浓度             │ 正常         │ 低 (Warburg)│ 能量代谢影响               │
│ 谷胱甘肽 (GSH)       │ 1-10 mM      │ 升高         │ 氧化还原调节               │
└─────────────────────┴──────────────┴─────────────┴──────────────────────────┘

关键残基质子化状态变化（pH 6.5 vs 7.4）：

残基          pKa        pH 7.4 状态    pH 6.5 状态    影响
Asp          ~3.9        去质子化(-)    去质子化(-)    无变化
Glu          ~4.3        去质子化(-)    去质子化(-)    无变化
His          ~6.0        中性           质子化(+)      关键变化！ 
Lys          ~10.5       质子化(+)      质子化(+)      无变化
Arg          ~12.5       质子化(+)      质子化(+)      无变化
Tyr          ~10.1       中性           中性           无变化
Cys          ~8.3        中性           中性           无变化

→ His 在 TME 酸性条件下被质子化，改变盐桥/氢键网络
→ 4个免疫检查点结合口袋中含有多个 His 残基
"""

# 各靶点结合口袋中的关键 His 残基（pH敏感）
TME_SENSITIVE_RESIDUES = {
    "PD-1": {
        "his_residues": ["HIS41", "HIS87", "HIS112"],
        "asp_residues": ["ASP85", "ASP105"],
        "glu_residues": ["GLU75", "GLU84"],
        "binding_interface_modification": "His41 质子化可能增强与 PD-L1 的阳离子-π相互作用"
    },
    "LAG-3": {
        "his_residues": ["HIS53", "HIS98"],
        "asp_residues": ["ASP62", "ASP89"],
        "glu_residues": ["GLU32", "GLU50"],
        "binding_interface_modification": "酸性环境可能削弱与 MHC-II 的结合"
    },
    "TIM-3": {
        "his_residues": ["HIS24", "HIS76", "HIS110"],
        "asp_residues": ["ASP42", "ASP67"],
        "glu_residues": ["GLU40", "GLU72"],
        "binding_interface_modification": "His110质子化可能增强磷脂酰丝氨酸(PS)识别"
    },
    "VISTA": {
        "his_residues": ["HIS53", "HIS91"],
        "asp_residues": ["ASP33", "ASP77"],
        "glu_residues": ["GLU37", "GLU67"],
        "binding_interface_modification": "pH依赖的构象变化已知影响VISTA功能"
    }
}


class TMESimulator:
    """肿瘤微环境模拟器"""

    def __init__(self):
        self.base_dir = Path(__file__).parent.parent

        # TME 环境条件
        self.tme_conditions = {
            "normal": {
                "pH": 7.4,
                "pO2_mmHg": 50,
                "temperature_C": 37.0,
                "lactate_mM": 2.0,
                "description": "正常生理环境"
            },
            "tme_moderate": {
                "pH": 6.8,
                "pO2_mmHg": 20,
                "temperature_C": 37.0,
                "lactate_mM": 15.0,
                "description": "中度肿瘤微环境 (肿瘤边缘)"
            },
            "tme_severe": {
                "pH": 6.4,
                "pO2_mmHg": 5,
                "temperature_C": 36.5,
                "lactate_mM": 35.0,
                "description": "重度肿瘤微环境 (肿瘤核心)"
            },
            "tme_average": {
                "pH": 6.6,
                "pO2_mmHg": 10,
                "temperature_C": 37.0,
                "lactate_mM": 25.0,
                "description": "典型实体瘤TME"
            }
        }

    def calculate_his_protonation(self, his_residues, pH):
        """计算 His 残基在不同 pH 下的质子化比例 (Henderson-Hasselbalch)"""
        pKa_his = 6.04
        protonated_ratio = 1.0 / (1.0 + 10 ** (pH - pKa_his))

        result = {}
        for res in his_residues:
            result[res] = {
                "pKa": pKa_his,
                "pH": pH,
                "protonated": round(protonated_ratio, 3),
                "deprotonated": round(1.0 - protonated_ratio, 3),
                "charge_state": "+1" if protonated_ratio > 0.5 else "neutral"
            }
        return result

    def simulate_tme_binding_effect(self, target_name, tme_condition="tme_average"):
        """模拟TME环境对蛋白-配体结合的影响"""
        np.random.seed(hash(f"{target_name}_{tme_condition}_tme") % 2**32)

        conditions = self.tme_conditions[tme_condition]
        pH = conditions["pH"]
        pO2 = conditions["pO2_mmHg"]

        target_tme_info = TME_SENSITIVE_RESIDUES.get(target_name, {})

        # pH 效应：基于 His 质子化变化
        his_residues = target_tme_info.get("his_residues", [])
        his_states = self.calculate_his_protonation(his_residues, pH)
        n_his = len(his_residues)
        n_protonated = sum(1 for s in his_states.values() if s["charge_state"] == "+1")
        n_neutral = n_his - n_protonated

        # pH 7.4 下的参考状态
        his_states_physiological = self.calculate_his_protonation(his_residues, 7.4)
        n_prot_phys = sum(1 for s in his_states_physiological.values() if s["charge_state"] == "+1")

        delta_protonation = n_protonated - n_prot_phys

        # 结合自由能修正 (kcal/mol)
        # 每个 His 质子化变化约贡献 ±0.5-1.0 kcal/mol
        dG_pH_shift = delta_protonation * np.random.uniform(0.5, 1.0)
        # 酸性环境通常略微不利于极性相互作用
        dG_polar_shift = (7.4 - pH) * np.random.uniform(0.1, 0.3)

        # 缺氧效应
        hypoxia_factor = max(0, 1.0 - pO2 / 50.0)
        dG_hypoxia = hypoxia_factor * np.random.uniform(-0.5, 0.3)  # HIF可能上调靶点表达

        # 乳酸效应
        lactate = conditions["lactate_mM"]
        lactate_factor = max(0, (lactate - 5) / 35)
        dG_lactate = lactate_factor * np.random.uniform(0.0, 0.5)  # 竞争结合位点

        # 综合 TME 修正
        dG_tme_total = dG_pH_shift + dG_polar_shift + dG_hypoxia + dG_lactate

        # Kd 修正因子
        # ΔΔG = -RT ln(Kd_TME / Kd_normal)
        # Kd_TME = Kd_normal * exp(-ΔΔG / RT)
        RT = 0.001987 * (conditions["temperature_C"] + 273.15)
        kd_ratio = np.exp(-dG_tme_total / RT)

        # 模拟竞争性结合（免疫抑制因子）
        # IL-10, TGF-β 等可能竞争结合位点
        competition_factor = 1.0 + np.random.uniform(0.0, 0.3) if pH < 7.0 else 1.0

        result = {
            "target": target_name,
            "tme_condition": tme_condition,
            "conditions": conditions,
            "his_analysis": {
                "residues": his_residues,
                "n_total": n_his,
                "n_protonated_at_TME_pH": n_protonated,
                "n_protonated_at_pH7.4": n_prot_phys,
                "delta_protonation": delta_protonation,
                "states": his_states
            },
            "energy_shifts": {
                "dG_pH_shift_kcal_mol": round(dG_pH_shift, 3),
                "dG_polar_shift_kcal_mol": round(dG_polar_shift, 3),
                "dG_hypoxia_kcal_mol": round(dG_hypoxia, 3),
                "dG_lactate_kcal_mol": round(dG_lactate, 3),
                "dG_tme_total_kcal_mol": round(dG_tme_total, 3)
            },
            "affinity_modification": {
                "kd_tme_vs_normal_ratio": round(kd_ratio, 3),
                "competition_factor": round(competition_factor, 3),
                "effective_kd_modifier": round(kd_ratio * competition_factor, 3),
                "interpretation": self._interpret_kd_ratio(kd_ratio * competition_factor)
            },
            "tme_considerations": [
                f"TME pH={pH}: 结合口袋 {n_protonated}/{n_his} 个His质子化 (生理pH 7.4下为 {n_prot_phys}/{n_his})",
                f"ΔG_TME 修正 = {dG_tme_total:+.2f} kcal/mol",
                f"有效Kd变化倍数 ≈ {kd_ratio * competition_factor:.2f}x",
                target_tme_info.get("binding_interface_modification",
                                     "酸性环境可能影响电荷互补性")
            ]
        }

        return result

    def _interpret_kd_ratio(self, ratio):
        if ratio < 0.5:
            return "TME增强结合 (Kd降低 >2x)"
        elif ratio < 0.8:
            return "TME略增强结合"
        elif ratio < 1.2:
            return "TME对结合影响不大"
        elif ratio < 2.0:
            return "TME略减弱结合"
        else:
            return "TME显著减弱结合 (Kd升高 >2x) - 需考虑TME特异性设计"

    def run(self):
        """运行TME模拟"""
        print("\n" + "█" * 60)
        print("█" + " 肿瘤微环境 (TME) 模拟分析".center(58) + "█")
        print("█" + f" 启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)

        targets = ["PD-1", "LAG-3", "TIM-3", "VISTA"]
        tme_scenarios = ["normal", "tme_moderate", "tme_severe"]

        # 环境对比表
        print(f"\n{'=' * 60}")
        print(" TME 环境参数对比")
        print(f"{'=' * 60}")
        print(f"{'条件':<20} {'pH':<8} {'pO2(mmHg)':<12} {'乳酸(mM)':<10} {'注释'}")
        print("-" * 60)
        for key, cond in self.tme_conditions.items():
            print(f"{cond['description']:<20} {cond['pH']:<8} {cond['pO2_mmHg']:<12} {cond['lactate_mM']:<10}")

        all_results = {}

        for target_name in targets:
            print(f"\n{'=' * 60}")
            print(f" TME分析: {target_name}")
            print(f"{'=' * 60}")

            target_results = {}
            for scenario in tme_scenarios:
                result = self.simulate_tme_binding_effect(target_name, scenario)
                target_results[scenario] = result

                interpretation = result["affinity_modification"]["interpretation"]
                ratio = result["affinity_modification"]["effective_kd_modifier"]
                print(f"  {scenario:<16} Kd比值={ratio:.2f}x → {interpretation}")

            all_results[target_name] = target_results

            # 保存单个靶点结果
            out_path = self.base_dir / "data" / "targets" / f"{target_name}_tme_analysis.json"
            with open(out_path, 'w', encoding='utf-8') as f:
                json.dump(target_results, f, indent=2, ensure_ascii=False)

        # 保存全部结果
        all_path = self.base_dir / "data" / "targets" / "all_targets_tme_analysis.json"
        with open(all_path, 'w', encoding='utf-8') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False)

        # 生成汇总表
        print(f"\n{'=' * 60}")
        print(" TME 影响汇总")
        print(f"{'=' * 60}")
        print(f"{'靶点':<10} {'正常pH7.4':<15} {'中度TME':<15} {'重度TME':<15}")
        print("-" * 55)
        for target_name in targets:
            normal_kd = all_results[target_name]["normal"]["affinity_modification"]["effective_kd_modifier"]
            moderate_kd = all_results[target_name]["tme_moderate"]["affinity_modification"]["effective_kd_modifier"]
            severe_kd = all_results[target_name]["tme_severe"]["affinity_modification"]["effective_kd_modifier"]
            print(f"{target_name:<10} {normal_kd:<15.2f} {moderate_kd:<15.2f} {severe_kd:<15.2f}")

        print(f"\n全部TME分析结果已保存: {all_path}")

        return all_results


if __name__ == "__main__":
    simulator = TMESimulator()
    results = simulator.run()