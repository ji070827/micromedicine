#!/usr/bin/env python3
"""
多免疫检查点小分子AI筛选系统 - Web可视化仪表盘
Flask + Chart.js 交互式数据可视化展示
"""

import os
import sys
import json
import threading
import numpy as np
import yaml
from pathlib import Path
from datetime import datetime
from flask import Flask, render_template, jsonify, request, send_file
from collections import defaultdict

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

import warnings
warnings.filterwarnings('ignore')

app = Flask(__name__, 
            template_folder='templates',
            static_folder='static')

# 全局状态
pipeline_status = {
    'running': False,
    'current_step': '',
    'progress': 0,
    'results': {},
    'route': None,
    'start_time': None
}

# 加载配置
config_path = Path(__file__).parent / "config" / "config.yaml"
with open(config_path, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)


def run_full_pipeline_route1():
    """运行路线一的完整管线"""
    global pipeline_status
    pipeline_status['running'] = True
    pipeline_status['route'] = 1
    pipeline_status['start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    try:
        # Step 1: 数据预处理
        pipeline_status['current_step'] = '数据预处理'
        pipeline_status['progress'] = 10
        from scripts.data_preprocess import DataPreprocessor
        preprocessor = DataPreprocessor()
        preprocess_results = preprocessor.run()
        pipeline_status['results']['preprocess'] = {
            'protein_targets': list(preprocess_results['protein_results'].keys()),
            'compound_count': len(preprocess_results['compounds_df']),
            'physchem_stats': preprocess_results['physchem_stats'],
            'viz_data': preprocess_results['viz_data']
        }
        
        # Step 2: DiffDock 批量对接
        pipeline_status['current_step'] = 'DiffDock 批量对接初筛'
        pipeline_status['progress'] = 25
        from scripts.diffdock_batch_run import DiffDockBatchRunner
        dd_runner = DiffDockBatchRunner()
        dd_results = dd_runner.run()
        
        # 加载可视化数据
        dd_viz_path = Path(__file__).parent / "results" / "diffdock" / "visualization_data.json"
        dd_viz = json.load(open(dd_viz_path, 'r', encoding='utf-8')) if dd_viz_path.exists() else {}
        pipeline_status['results']['diffdock'] = dd_viz
        
        # Step 3: 初筛过滤
        pipeline_status['current_step'] = '初筛排序与过滤'
        pipeline_status['progress'] = 45
        from scripts.primary_screen_filter import PrimaryScreenFilter
        filter = PrimaryScreenFilter()
        filter_results = filter.run()
        
        screen_viz_path = Path(__file__).parent / "results" / "primary_screen" / "visualization_data.json"
        screen_viz = json.load(open(screen_viz_path, 'r', encoding='utf-8')) if screen_viz_path.exists() else {}
        pipeline_status['results']['primary_screen'] = screen_viz
        
        # Step 4: AlphaFold3 复合物预测
        pipeline_status['current_step'] = 'AlphaFold3 复合物精细模拟'
        pipeline_status['progress'] = 65
        from scripts.af3_complex_prediction import AlphaFold3Runner
        af3_runner = AlphaFold3Runner()
        af3_results = af3_runner.run()
        
        af3_viz_path = Path(__file__).parent / "results" / "alphafold3" / "visualization_data.json"
        af3_viz = json.load(open(af3_viz_path, 'r', encoding='utf-8')) if af3_viz_path.exists() else {}
        pipeline_status['results']['alphafold3'] = af3_viz
        
        # Step 5: 相互作用分析
        pipeline_status['current_step'] = '蛋白-配体相互作用分析'
        pipeline_status['progress'] = 70
        from scripts.interaction_analysis import InteractionAnalyzer
        analyzer = InteractionAnalyzer()
        ia_results = analyzer.run()
        
        ia_viz_path = Path(__file__).parent / "results" / "alphafold3" / "interaction_viz_data.json"
        ia_viz = json.load(open(ia_viz_path, 'r', encoding='utf-8')) if ia_viz_path.exists() else {}
        pipeline_status['results']['interaction'] = ia_viz
        
        # Step 5.5: 功能预测（竞争性结合 + 选择性 + ADME/Tox）
        pipeline_status['current_step'] = '功能预测：竞争性结合+选择性+ADME/Tox'
        pipeline_status['progress'] = 78
        from scripts.competitive_binding import CompetitiveBindingPredictor
        cb_predictor = CompetitiveBindingPredictor()
        cb_results = cb_predictor.run()
        
        from scripts.selectivity_analysis import SelectivityAnalyzer
        sel_analyzer = SelectivityAnalyzer()
        sel_results = sel_analyzer.run()
        
        from scripts.adme_predictor import ADMEPredictor
        adme_predictor = ADMEPredictor()
        adme_results = adme_predictor.run()
        
        # 加载功能预测可视化数据
        cb_viz_path = Path(__file__).parent / "results" / "alphafold3" / "competitive_binding_viz.json"
        cb_viz = json.load(open(cb_viz_path, 'r', encoding='utf-8')) if cb_viz_path.exists() else {}
        pipeline_status['results']['competitive_binding'] = cb_viz
        
        # Step 6: 终选排序（七维加权打分）
        pipeline_status['current_step'] = '终选排序与报告生成（七维加权）'
        pipeline_status['progress'] = 90
        from scripts.final_ranking import FinalRanking
        ranker = FinalRanking()
        ranking_results = ranker.run()
        
        final_viz_path = Path(__file__).parent / "results" / "final_report" / "visualization_data.json"
        final_viz = json.load(open(final_viz_path, 'r', encoding='utf-8')) if final_viz_path.exists() else {}
        pipeline_status['results']['final_ranking'] = final_viz
        
        # 加载最终报告
        report_path = Path(__file__).parent / "results" / "final_report" / "final_report.json"
        if report_path.exists():
            final_report = json.load(open(report_path, 'r', encoding='utf-8'))
            pipeline_status['results']['final_report'] = final_report
        
        pipeline_status['progress'] = 100
        pipeline_status['current_step'] = '完成'
        
    except Exception as e:
        pipeline_status['current_step'] = f'错误: {str(e)}'
        pipeline_status['error'] = str(e)
        print(f"Pipeline error: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        pipeline_status['running'] = False


def run_full_pipeline_route2():
    """运行路线二的完整管线"""
    global pipeline_status
    pipeline_status['running'] = True
    pipeline_status['route'] = 2
    pipeline_status['start_time'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    try:
        # Step 1: 活性模型训练
        pipeline_status['current_step'] = '活性特征提取与模型训练'
        pipeline_status['progress'] = 15
        from scripts.activity_model_train import ActivityModelTrainer
        trainer = ActivityModelTrainer()
        model_results = trainer.run()
        
        model_viz_path = Path(__file__).parent / "data" / "activity_dataset" / "model_viz_data.json"
        model_viz = json.load(open(model_viz_path, 'r', encoding='utf-8')) if model_viz_path.exists() else {}
        pipeline_status['results']['activity_model'] = model_viz
        
        # Step 2: 全新分子生成（TargetDiff 口袋感知生成，降级到 RDKit 组合化学）
        pipeline_status['current_step'] = '全新分子生成与初筛 (TargetDiff)'
        pipeline_status['progress'] = 45
        from scripts.targetdiff_generate import TargetDiffGenerator
        generator = TargetDiffGenerator()
        gen_results = generator.run()
        
        gen_viz_path = Path(__file__).parent / "data" / "activity_dataset" / "generation_viz_data.json"
        gen_viz = json.load(open(gen_viz_path, 'r', encoding='utf-8')) if gen_viz_path.exists() else {}
        pipeline_status['results']['molecule_generation'] = gen_viz
        
        # Step 3: 数据预处理
        pipeline_status['current_step'] = '数据预处理'
        pipeline_status['progress'] = 55
        from scripts.data_preprocess import DataPreprocessor
        preprocessor = DataPreprocessor()
        preprocess_results = preprocessor.run()
        pipeline_status['results']['preprocess'] = {
            'protein_targets': list(preprocess_results['protein_results'].keys()),
            'compound_count': len(preprocess_results['compounds_df']),
            'physchem_stats': preprocess_results['physchem_stats'],
            'viz_data': preprocess_results['viz_data']
        }
        
        # Step 4: DiffDock对接 + 后续流程
        pipeline_status['current_step'] = 'DiffDock 对接与后续筛选'
        pipeline_status['progress'] = 70
        from scripts.diffdock_batch_run import DiffDockBatchRunner
        dd_runner = DiffDockBatchRunner()
        dd_runner.run()
        
        from scripts.primary_screen_filter import PrimaryScreenFilter
        filter = PrimaryScreenFilter()
        filter.run()
        
        from scripts.af3_complex_prediction import AlphaFold3Runner
        af3_runner = AlphaFold3Runner()
        af3_runner.run()
        
        from scripts.interaction_analysis import InteractionAnalyzer
        analyzer = InteractionAnalyzer()
        analyzer.run()
        
        pipeline_status['current_step'] = '终选排序与报告生成'
        pipeline_status['progress'] = 90
        from scripts.final_ranking import FinalRanking
        ranker = FinalRanking()
        ranking_results = ranker.run()
        
        final_viz_path = Path(__file__).parent / "results" / "final_report" / "visualization_data.json"
        final_viz = json.load(open(final_viz_path, 'r', encoding='utf-8')) if final_viz_path.exists() else {}
        pipeline_status['results']['final_ranking'] = final_viz
        
        report_path = Path(__file__).parent / "results" / "final_report" / "final_report.json"
        if report_path.exists():
            final_report = json.load(open(report_path, 'r', encoding='utf-8'))
            pipeline_status['results']['final_report'] = final_report
        
        pipeline_status['progress'] = 100
        pipeline_status['current_step'] = '完成'
        
    except Exception as e:
        pipeline_status['current_step'] = f'错误: {str(e)}'
        pipeline_status['error'] = str(e)
        print(f"Pipeline error: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        pipeline_status['running'] = False


def generate_dashboard_data():
    """生成仪表盘所需的聚合数据"""
    targets = ['PD-1', 'LAG-3', 'TIM-3', 'VISTA']
    target_colors = config['visualization']['target_colors']
    
    data = {
        'targets': targets,
        'target_colors': target_colors,
        'target_info': {
            t: {
                'name': t,
                'description': config['targets'][t]['description'],
                'pdb_id': config['targets'][t]['pdb_id'],
                'functional_domain': config['targets'][t]['functional_domain'],
                'n_binding_residues': len(config['targets'][t]['binding_site_residues'])
            }
            for t in targets
        },
        'weights': config['screening']['final_ranking']['weights'],
        'screening_params': {
            'diffdock': config['screening']['diffdock'],
            'primary_filter_top_n': config['screening']['primary_filter']['top_n_per_target'],
            'alphafold3': config['screening']['alphafold3']
        },
        'route2_params': config['route2']
    }
    
    # 尝试加载已保存的各类可视化数据
    # DiffDock数据
    dd_viz_path = Path(__file__).parent / "results" / "diffdock" / "visualization_data.json"
    if dd_viz_path.exists():
        data['diffdock'] = json.load(open(dd_viz_path, 'r', encoding='utf-8'))
        data['has_diffdock'] = True
    else:
        data['has_diffdock'] = False
    
    # 初筛数据
    screen_viz_path = Path(__file__).parent / "results" / "primary_screen" / "visualization_data.json"
    if screen_viz_path.exists():
        data['primary_screen'] = json.load(open(screen_viz_path, 'r', encoding='utf-8'))
        data['has_primary_screen'] = True
    else:
        data['has_primary_screen'] = False
    
    # AlphaFold3数据
    af3_viz_path = Path(__file__).parent / "results" / "alphafold3" / "visualization_data.json"
    if af3_viz_path.exists():
        data['alphafold3'] = json.load(open(af3_viz_path, 'r', encoding='utf-8'))
        data['has_alphafold3'] = True
    else:
        data['has_alphafold3'] = False
    
    # 相互作用数据
    ia_viz_path = Path(__file__).parent / "results" / "alphafold3" / "interaction_viz_data.json"
    if ia_viz_path.exists():
        data['interaction'] = json.load(open(ia_viz_path, 'r', encoding='utf-8'))
        data['has_interaction'] = True
    else:
        data['has_interaction'] = False
    
    # 终选数据
    final_viz_path = Path(__file__).parent / "results" / "final_report" / "visualization_data.json"
    if final_viz_path.exists():
        data['final_ranking'] = json.load(open(final_viz_path, 'r', encoding='utf-8'))
        data['has_final_ranking'] = True
    else:
        data['has_final_ranking'] = False
    
    # 最终报告
    report_path = Path(__file__).parent / "results" / "final_report" / "final_report.json"
    if report_path.exists():
        data['final_report'] = json.load(open(report_path, 'r', encoding='utf-8'))
        data['has_final_report'] = True
    else:
        data['has_final_report'] = False
    
    # 预处理数据
    lib_viz_path = Path(__file__).parent / "data" / "library" / "visualization_data.json"
    if lib_viz_path.exists():
        data['preprocess'] = json.load(open(lib_viz_path, 'r', encoding='utf-8'))
        data['has_preprocess'] = True
    else:
        data['has_preprocess'] = False
    
    # 活性模型数据（路线二）
    model_viz_path = Path(__file__).parent / "data" / "activity_dataset" / "model_viz_data.json"
    if model_viz_path.exists():
        data['activity_model'] = json.load(open(model_viz_path, 'r', encoding='utf-8'))
        data['has_activity_model'] = True
    else:
        data['has_activity_model'] = False
    
    # 分子生成数据（路线二）
    gen_viz_path = Path(__file__).parent / "data" / "activity_dataset" / "generation_viz_data.json"
    if gen_viz_path.exists():
        data['molecule_generation'] = json.load(open(gen_viz_path, 'r', encoding='utf-8'))
        data['has_molecule_generation'] = True
    else:
        data['has_molecule_generation'] = False
    
    return data


# ==================== 路由 ====================

@app.route('/')
def index():
    """主页 - 控制台仪表盘"""
    return render_template('index.html')


@app.route('/route1')
def route1_page():
    """路线一：已知库筛选"""
    return render_template('route1.html')


@app.route('/route2')
def route2_page():
    """路线二：全新分子设计"""
    return render_template('route2.html')


@app.route('/results')
def results_page():
    """结果总览"""
    return render_template('results.html')


@app.route('/api/status')
def get_status():
    """获取管线运行状态"""
    return jsonify(pipeline_status)


@app.route('/api/dashboard_data')
def get_dashboard_data():
    """获取仪表盘数据"""
    data = generate_dashboard_data()
    return jsonify(data)


@app.route('/api/run_pipeline', methods=['POST'])
def run_pipeline():
    """运行完整筛选管线"""
    global pipeline_status
    
    if pipeline_status['running']:
        return jsonify({'error': 'Pipeline is already running'}), 400
    
    route = request.json.get('route', 1)
    
    if route == 1:
        thread = threading.Thread(target=run_full_pipeline_route1, daemon=True)
    else:
        thread = threading.Thread(target=run_full_pipeline_route2, daemon=True)
    
    thread.start()
    
    return jsonify({'message': f'Pipeline Route {route} started', 'route': route})


@app.route('/api/run_step', methods=['POST'])
def run_single_step():
    """运行单个步骤"""
    step = request.json.get('step', '')
    try:
        result = {}
        if step == 'preprocess':
            from scripts.data_preprocess import DataPreprocessor
            p = DataPreprocessor()
            r = p.run()
            result = {'status': 'ok', 'compounds': len(r['compounds_df']), 'targets': list(r['protein_results'].keys())}
        elif step == 'diffdock':
            from scripts.diffdock_batch_run import DiffDockBatchRunner
            r = DiffDockBatchRunner()
            r.run()
            result = {'status': 'ok'}
        elif step == 'primary_screen':
            from scripts.primary_screen_filter import PrimaryScreenFilter
            r = PrimaryScreenFilter()
            r.run()
            result = {'status': 'ok'}
        elif step == 'alphafold3':
            from scripts.af3_complex_prediction import AlphaFold3Runner
            r = AlphaFold3Runner()
            r.run()
            result = {'status': 'ok'}
        elif step == 'interaction':
            from scripts.interaction_analysis import InteractionAnalyzer
            r = InteractionAnalyzer()
            r.run()
            result = {'status': 'ok'}
        elif step == 'final_ranking':
            from scripts.final_ranking import FinalRanking
            r = FinalRanking()
            r.run()
            result = {'status': 'ok'}
        elif step == 'activity_model':
            from scripts.activity_model_train import ActivityModelTrainer
            r = ActivityModelTrainer()
            r.run()
            result = {'status': 'ok'}
        elif step == 'molecule_generation':
            from scripts.targetdiff_generate import TargetDiffGenerator
            r = TargetDiffGenerator()
            r.run()
            result = {'status': 'ok'}
        else:
            result = {'status': 'error', 'message': f'Unknown step: {step}'}
        
        return jsonify(result)
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)}), 500


@app.route('/api/scores/<target>')
def get_target_scores(target):
    """获取特定靶点的详细得分数据"""
    # 模拟各维度得分数据
    np.random.seed(hash(target) % 2**32)
    
    n_candidates = 75
    scores = []
    for i in range(n_candidates):
        scores.append({
            'rank': i + 1,
            'mol_id': f'CMPD_{i+1:06d}',
            'docking_score': round(np.clip(np.random.beta(3, 3), 0.1, 0.99), 3),
            'structure_confidence': round(np.clip(np.random.beta(5, 2), 0.2, 0.99), 3),
            'interaction_strength': round(np.clip(np.random.beta(4, 3), 0.1, 0.99), 3),
            'drug_likeness': round(np.random.beta(5, 2), 3),
            'binding_site_match': round(np.random.beta(4, 2.5), 3),
            'final_score': round(np.random.beta(5, 3) * 0.85 + 0.1, 3)
        })
    
    # 按final_score排序
    scores.sort(key=lambda x: x['final_score'], reverse=True)
    for i, s in enumerate(scores):
        s['rank'] = i + 1
        if s['final_score'] >= 0.75:
            s['priority'] = 'A'
        elif s['final_score'] >= 0.6:
            s['priority'] = 'B'
        elif s['final_score'] >= 0.45:
            s['priority'] = 'C'
        else:
            s['priority'] = 'D'
    
    return jsonify(scores)


@app.route('/api/3d_complex/<target>')
def get_3d_complex(target):
    """获取指定靶点的3D复合物结构数据（含所有候选药物，用于3Dmol.js可视化）"""
    if target not in config['targets']:
        return jsonify({'error': f'Unknown target: {target}'}), 404
    
    # 尝试加载已生成的3D结构
    complex_path = Path(__file__).parent / "data" / "targets" / f"{target}_3d_complex.json"
    if complex_path.exists():
        with open(complex_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # 检查是否是新格式（含drugs列表）
        if 'drugs' in data:
            return jsonify(data)
    
    # 如果没有预先计算或格式过旧，实时生成
    from scripts.generate_3d_complex import Complex3DGenerator
    gen = Complex3DGenerator()
    result = gen.generate_complex_with_drugs(target)
    if result:
        return jsonify(result)
    
    return jsonify({'error': 'Failed to generate complex structure'}), 500


@app.route('/api/3d_complex/<target>/<drug_id>')
def get_3d_drug_complex(target, drug_id):
    """获取指定靶点+指定药物的复合物PDB数据"""
    if target not in config['targets']:
        return jsonify({'error': f'Unknown target: {target}'}), 404
    
    # 先尝试从完整JSON中提取
    complex_path = Path(__file__).parent / "data" / "targets" / f"{target}_3d_complex.json"
    if complex_path.exists():
        with open(complex_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        if 'drugs' in data:
            for drug in data['drugs']:
                if drug['mol_id'] == drug_id:
                    return jsonify({
                        'target': target,
                        'mol_id': drug_id,
                        'pdb_string': drug['pdb_string'],
                        'name': drug['name'],
                        'mol_type': drug.get('mol_type', ''),
                        'binding_affinity': drug.get('binding_affinity', 0),
                        'predicted_interactions': drug.get('predicted_interactions', {}),
                        'n_ligand_atoms': drug.get('n_ligand_atoms', 0),
                        'binding_pocket_center': data.get('binding_pocket_center', [0,0,0]),
                        'binding_pocket_residues': data.get('binding_pocket_residues', []),
                        'n_protein_residues': data.get('n_protein_residues', 0),
                    })
    
    # 如果文件不存在，实时生成
    from scripts.generate_3d_complex import Complex3DGenerator
    gen = Complex3DGenerator()
    result = gen.generate_complex_with_drugs(target)
    if result and 'drugs' in result:
        for drug in result['drugs']:
            if drug['mol_id'] == drug_id:
                return jsonify({
                    'target': target,
                    'mol_id': drug_id,
                    'pdb_string': drug['pdb_string'],
                    'name': drug['name'],
                    'mol_type': drug.get('mol_type', ''),
                    'binding_affinity': drug.get('binding_affinity', 0),
                    'predicted_interactions': drug.get('predicted_interactions', {}),
                    'n_ligand_atoms': drug.get('n_ligand_atoms', 0),
                    'binding_pocket_center': result.get('binding_pocket_center', [0,0,0]),
                    'binding_pocket_residues': result.get('binding_pocket_residues', []),
                    'n_protein_residues': result.get('n_protein_residues', 0),
                })
    
    return jsonify({'error': f'Drug {drug_id} not found for target {target}'}), 404


@app.route('/api/drug_list/<target>')
def get_drug_list(target):
    """获取指定靶点的所有候选药物列表（轻量）"""
    if target not in config['targets']:
        return jsonify({'error': f'Unknown target: {target}'}), 404
    
    complex_path = Path(__file__).parent / "data" / "targets" / f"{target}_3d_complex.json"
    if complex_path.exists():
        with open(complex_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        if 'drugs' in data:
            drugs_light = []
            for d in data['drugs']:
                drugs_light.append({
                    'mol_id': d['mol_id'],
                    'name': d['name'],
                    'mol_type': d.get('mol_type', ''),
                    'smiles_like': d.get('smiles_like', ''),
                    'binding_affinity': d.get('binding_affinity', 0),
                    'mw': d.get('mw', 0),
                    'drug_likeness': d.get('drug_likeness', 0),
                    'predicted_interactions': d.get('predicted_interactions', {}),
                    'n_ligand_atoms': d.get('n_ligand_atoms', 0),
                })
            return jsonify({'target': target, 'n_drugs': len(drugs_light), 'drugs': drugs_light})
    
    # 生成
    from scripts.generate_3d_complex import Complex3DGenerator
    gen = Complex3DGenerator()
    result = gen.generate_complex_with_drugs(target)
    if result and 'drugs' in result:
        drugs_light = []
        for d in result['drugs']:
            drugs_light.append({
                'mol_id': d['mol_id'],
                'name': d['name'],
                'mol_type': d.get('mol_type', ''),
                'smiles_like': d.get('smiles_like', ''),
                'binding_affinity': d.get('binding_affinity', 0),
                'mw': d.get('mw', 0),
                'drug_likeness': d.get('drug_likeness', 0),
                'predicted_interactions': d.get('predicted_interactions', {}),
                'n_ligand_atoms': d.get('n_ligand_atoms', 0),
            })
        return jsonify({'target': target, 'n_drugs': len(drugs_light), 'drugs': drugs_light})
    
    return jsonify({'error': f'Failed to get drug list for {target}'}), 500


@app.route('/api/all_drug_lists')
def get_all_drug_lists():
    """获取所有靶点的候选药物列表摘要"""
    targets = ['PD-1', 'LAG-3', 'TIM-3', 'VISTA']
    all_lists = {}
    
    drug_list_path = Path(__file__).parent / "data" / "targets" / "all_drug_candidates.json"
    if drug_list_path.exists():
        with open(drug_list_path, 'r', encoding='utf-8') as f:
            return jsonify(json.load(f))
    
    # 如果没有缓存文件，逐个生成
    from scripts.generate_3d_complex import Complex3DGenerator
    gen = Complex3DGenerator()
    
    for target in targets:
        complex_path = Path(__file__).parent / "data" / "targets" / f"{target}_3d_complex.json"
        if complex_path.exists():
            with open(complex_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            if 'drugs' in data:
                all_lists[target] = [
                    {
                        'mol_id': d['mol_id'],
                        'name': d['name'],
                        'mol_type': d.get('mol_type', ''),
                        'binding_affinity': d.get('binding_affinity', 0),
                        'mw': d.get('mw', 0),
                        'drug_likeness': d.get('drug_likeness', 0),
                        'predicted_interactions': d.get('predicted_interactions', {}),
                        'n_ligand_atoms': d.get('n_ligand_atoms', 0),
                    }
                    for d in data['drugs']
                ]
                continue
        
        result = gen.generate_complex_with_drugs(target)
        if result and 'drugs' in result:
            all_lists[target] = [
                {
                    'mol_id': d['mol_id'],
                    'name': d['name'],
                    'mol_type': d.get('mol_type', ''),
                    'binding_affinity': d.get('binding_affinity', 0),
                    'mw': d.get('mw', 0),
                    'drug_likeness': d.get('drug_likeness', 0),
                    'predicted_interactions': d.get('predicted_interactions', {}),
                    'n_ligand_atoms': d.get('n_ligand_atoms', 0),
                }
                for d in result['drugs']
            ]
    
    return jsonify(all_lists)


@app.route('/structure')
def structure_viewer():
    """3D结构查看器页面"""
    return render_template('structure.html')


@app.route('/api/upload_dataset', methods=['POST'])
def upload_dataset():
    """上传新的小分子数据集文件 (CSV/SMI/SDF)"""
    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
    
    # 验证扩展名
    ext = os.path.splitext(file.filename)[1].lower()
    if ext not in ['.csv', '.smi', '.sdf', '.txt']:
        return jsonify({'error': f'Unsupported format: {ext}. Supported: .csv, .smi, .sdf, .txt'}), 400
    
    # 保存到 upload_dir
    upload_dir = Path(__file__).parent / "data" / "library" / "user_uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)
    
    save_path = upload_dir / file.filename
    file.save(str(save_path))
    
    # 解析文件，提取基本信息
    n_compounds = 0
    format_type = ext.lstrip('.')
    
    try:
        if ext == '.csv':
            import pandas as pd
            df = pd.read_csv(save_path, nrows=5)
            n_compounds = sum(1 for _ in open(save_path, 'r', encoding='utf-8')) - 1  # rough count
            columns = list(df.columns)
        elif ext in ['.smi', '.txt']:
            with open(save_path, 'r', encoding='utf-8') as f:
                lines = [l.strip() for l in f if l.strip()]
                n_compounds = len(lines)
            columns = ['smiles']
        elif ext == '.sdf':
            with open(save_path, 'r', encoding='utf-8') as f:
                content = f.read()
                n_compounds = content.count('$$$$')
            columns = ['MOL_block']
    except Exception as e:
        return jsonify({'error': f'Failed to parse file: {str(e)}'}), 400
    
    # 确定数据集名称
    dataset_name = request.form.get('name', file.filename)
    
    result = {
        'status': 'ok',
        'filename': file.filename,
        'name': dataset_name,
        'format': format_type,
        'n_compounds': n_compounds,
        'columns': columns[:10],  # first 10 columns only
        'path': str(save_path),
    }
    
    # 自动切换为激活数据集
    config['data']['active_library'] = file.filename
    # 也保存到 library 目录一份（如果格式支持直接使用）
    lib_copy = Path(__file__).parent / "data" / "library" / file.filename
    if save_path != lib_copy:
        import shutil
        shutil.copy2(str(save_path), str(lib_copy))
        result['activated'] = True
        result['active_path'] = str(lib_copy)
    
    return jsonify(result)


@app.route('/api/list_datasets')
def list_datasets():
    """列出所有可用的小分子数据集"""
    library_dir = Path(__file__).parent / "data" / "library"
    upload_dir = library_dir / "user_uploads"
    
    datasets = []
    
    # 扫描 library 目录
    for pattern in ['*.csv', '*.smi', '*.sdf']:
        for fpath in library_dir.glob(pattern):
            if fpath.parent == upload_dir:
                continue  # 跳过 user_uploads 中的，会单独列出
            stat = fpath.stat()
            datasets.append({
                'filename': fpath.name,
                'format': fpath.suffix.lstrip('.'),
                'size_kb': round(stat.st_size / 1024, 1),
                'modified': datetime.fromtimestamp(stat.st_mtime).strftime('%Y-%m-%d %H:%M'),
                'source': 'library',
                'active': config.get('data', {}).get('active_library', '') == fpath.name,
            })
    
    # 扫描 user_uploads 目录
    if upload_dir.exists():
        for fpath in upload_dir.glob('*'):
            if fpath.suffix.lower() in ['.csv', '.smi', '.sdf', '.txt']:
                stat = fpath.stat()
                datasets.append({
                    'filename': fpath.name,
                    'format': fpath.suffix.lstrip('.'),
                    'size_kb': round(stat.st_size / 1024, 1),
                    'modified': datetime.fromtimestamp(stat.st_mtime).strftime('%Y-%m-%d %H:%M'),
                    'source': 'user_uploads',
                    'active': config.get('data', {}).get('active_library', '') == fpath.name,
                })
    
    # 按修改时间倒序
    datasets.sort(key=lambda x: x['modified'], reverse=True)
    
    active_lib = config.get('data', {}).get('active_library', 'pubchem_all_targets.csv')
    return jsonify({
        'active_library': active_lib,
        'n_datasets': len(datasets),
        'datasets': datasets,
    })


@app.route('/api/switch_dataset', methods=['POST'])
def switch_dataset():
    """切换当前激活的数据集"""
    data = request.get_json() or {}
    filename = data.get('filename', '')
    
    if not filename:
        return jsonify({'error': 'filename is required'}), 400
    
    # 检查文件是否存在
    library_dir = Path(__file__).parent / "data" / "library"
    upload_dir = library_dir / "user_uploads"
    
    target_path = library_dir / filename
    if not target_path.exists():
        target_path = upload_dir / filename
    if not target_path.exists():
        return jsonify({'error': f'Dataset not found: {filename}'}), 404
    
    # 更新配置
    config['data']['active_library'] = filename
    
    # 如果文件只在 user_uploads 中，复制到 library
    if not (library_dir / filename).exists() and (upload_dir / filename).exists():
        import shutil
        shutil.copy2(str(upload_dir / filename), str(library_dir / filename))
    
    return jsonify({
        'status': 'ok',
        'active_library': filename,
        'message': f'Switched to dataset: {filename}',
    })


@app.route('/api/config')
def get_config():
    """获取当前配置"""
    return jsonify(config)


if __name__ == '__main__':
    print("\n" + "=" * 60)
    print("  多免疫检查点小分子AI筛选系统")
    print("  Immuno-Checkpoint Small Molecule AI Screening")
    print("=" * 60)
    print("\n  正在启动 Web 仪表盘...")
    print("  本机访问: http://127.0.0.1:5050")
    print("  局域网访问: http://<本机IP>:5050")
    print("  按 Ctrl+C 停止服务器")
    print("\n" + "=" * 60 + "\n")
    # 注意: debug=False 且 use_reloader=False 避免 Windows 下
    # 多进程/端口占用问题，确保服务器稳定可连接。
    app.run(host='0.0.0.0', port=5050, debug=False, use_reloader=False)
