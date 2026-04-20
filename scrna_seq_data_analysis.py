#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRNA-seq Data Analysis Pipeline
Containerized version for reproducible single-cell analysis
"""

import os
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import matplotlib.pyplot as plt

# Non-interactive backend для сохранения графиков в файл
plt.switch_backend('Agg')


DATA_URL = os.getenv('DATA_URL', 'https://datasets.cellxgene.cziscience.com/b9cbe943-ad26-4cac-8798-6453b80834bf.h5ad')
DATA_FILE = os.getenv('DATA_FILE', 'HumanBrain_NuclAccumb.h5ad')
DATA_DIR = os.getenv('DATA_DIR', '/data')
OUTPUT_DIR = os.getenv('OUTPUT_DIR', '/output')
MIN_COUNTS_PER_GENE = int(os.getenv('MIN_COUNTS_PER_GENE', '1000'))
MIN_COUNTS_PER_CELL = int(os.getenv('MIN_COUNTS_PER_CELL', '15000'))
APP_PORT = int(os.getenv('APP_PORT', '8080'))
WEB_MODE = os.getenv('WEB_MODE', 'false').lower() == 'true'


def setup_directories():
    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def download_data():
    data_path = Path(DATA_DIR) / DATA_FILE
    if not data_path.exists():
        print(f"[INFO] Downloading data from {DATA_URL}...")
        import urllib.request
        urllib.request.urlretrieve(DATA_URL, str(data_path))
        print(f"[OK] Data saved: {data_path}")
    else:
        print(f"[INFO] Data already exists: {data_path}")
    return str(data_path)

def load_data(filepath):
    print(f"[INFO] Loading AnnData from {filepath}...")
    adata = sc.read_h5ad(filepath)
    print(f"[OK] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata

def run_qc(adata):
    print("[INFO] Running QC metrics...")
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # Убираем гены с аномально высокими счетами
    adata = adata[:, adata.var["total_counts"] < 1e7]
    return adata

def remove_doublets(adata):
    print("[INFO] Running Scrublet for doublet detection...")
    try:
        sce.pp.scrublet(adata)
        adata = adata[~adata.obs['predicted_doublet'], :]
        print(f"[OK] After doublet removal: {adata.n_obs} cells")
    except Exception as e:
        print(f"[WARN] Scrublet skipped: {e}")
    return adata

def filter_low_counts(adata, min_gene, min_cell):
    print(f"[INFO] Filtering: genes≥{min_gene}, cells≥{min_cell}")
    cell_counts = adata.obs['total_UMIs'] if 'total_UMIs' in adata.obs else adata.obs['n_counts']
    gene_counts = adata.var['total_counts']
    adata = adata[:, gene_counts >= min_gene]
    adata = adata[cell_counts >= min_cell, :]
    print(f"[OK] After filtering: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata

def run_dimensionality_reduction(adata):
    print("[INFO] Running PCA...")
    sc.tl.pca(adata)
    print("[INFO] Running UMAP...")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata

def save_outputs(adata, output_dir):
    print(f"[INFO] Saving results to {output_dir}...")
    # Нормализация
    sc.pp.log1p(adata)
    # Сохранение обработанных данных
    processed_path = Path(output_dir) / "processed_adata.h5ad"
    adata.write(processed_path)
    print(f"[OK] Saved: {processed_path}")
    
    # Сохранение графиков
    color_cols = [c for c in ['supercluster_term', 'donor_id', 'cell_type'] if c in adata.obs.columns]
    for col in color_cols[:2]:
        try:
            sc.pl.umap(adata, color=col, save=f"_{col}.png", show=False)
        except:
            pass
    try:
        sc.pl.pca_variance_ratio(adata, n_pcs=20, save="_pca.png", show=False)
    except:
        pass
    print(f"[OK] Plots saved to {output_dir}")

def run_web_server(output_dir, port):
    from flask import Flask, jsonify, send_from_directory
    app = Flask(__name__)
    
    @app.route('/')
    def index():
        return jsonify({
            "service": "scRNA-seq Analysis Results",
            "status": "ready",
            "output_dir": output_dir
        })
    
    @app.route('/files')
    def list_files():
        files = [f.name for f in Path(output_dir).glob('*') if f.is_file()]
        return jsonify({"files": files})
    
    @app.route('/plot/<path:filename>')
    def serve_plot(filename):
        return send_from_directory(output_dir, filename)
    
    print(f"[INFO] Starting web server on 0.0.0.0:{port}")
    app.run(host='0.0.0.0', port=port)

def main():
    print("scRNA-seq Analysis Pipeline")
    print(f"Config: DATA_DIR={DATA_DIR}, OUTPUT_DIR={OUTPUT_DIR}, WEB_MODE={WEB_MODE}")
    
    setup_directories()
    data_path = download_data()
    adata = load_data(data_path)
    
    # Analysis pipeline
    adata = run_qc(adata)
    adata = remove_doublets(adata)
    adata = filter_low_counts(adata, MIN_COUNTS_PER_GENE, MIN_COUNTS_PER_CELL)
    adata = run_dimensionality_reduction(adata)
    save_outputs(adata, OUTPUT_DIR)
    
    print("Analysis complete")
    
    if WEB_MODE:
        run_web_server(OUTPUT_DIR, APP_PORT)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
