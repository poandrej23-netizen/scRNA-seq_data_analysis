#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scRNA-seq Data Analysis Pipeline
Containerized version with interactive output & demo mode
"""

import os
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import matplotlib.pyplot as plt


plt.switch_backend('Agg')


DATA_URL = os.getenv('DATA_URL', 'https://datasets.cellxgene.cziscience.com/b9cbe943-ad26-4cac-8798-6453b80834bf.h5ad')
DATA_FILE = os.getenv('DATA_FILE', 'HumanBrain_NuclAccumb.h5ad')
DATA_DIR = os.getenv('DATA_DIR', '/data')
OUTPUT_DIR = os.getenv('OUTPUT_DIR', '/output')
MIN_COUNTS_PER_GENE = int(os.getenv('MIN_COUNTS_PER_GENE', '1000'))
MIN_COUNTS_PER_CELL = int(os.getenv('MIN_COUNTS_PER_CELL', '15000'))
APP_PORT = int(os.getenv('APP_PORT', '8080'))
WEB_MODE = os.getenv('WEB_MODE', 'false').lower() == 'true'
DEMO_MODE = os.getenv('DEMO_MODE', 'false').lower() == 'true'
DEMO_CELLS = int(os.getenv('DEMO_CELLS', '2000'))  # subsample для демо


def print_step(msg: str, demo: bool = False):
    prefix = "[DEMO] " if demo else "."
    print(f"{prefix}{msg}", flush=True)

def download_with_progress(url: str, filepath: str, chunk_size: int = 8192):
    from urllib.error import URLError
    
    print(f"[INFO] Downloading: {url}", flush=True)
    
    req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urlopen(req, timeout=30) as response:
            total_size = int(response.headers.get('Content-Length', 0))
            downloaded = 0
            start_time = time.time()
            
            with open(filepath, 'wb') as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    
                    # Прогресс-бар в терминал
                    if total_size:
                        percent = downloaded / total_size * 100
                        elapsed = time.time() - start_time
                        speed = downloaded / elapsed / 1024 / 1024 if elapsed > 0 else 0
                        bar_len = 40
                        filled = int(bar_len * percent / 100)
                        bar = '+' * filled + '.' * (bar_len - filled)
                        print(f"[{bar}] {percent}% ({downloaded/1e6}/{total_size/1e6} MB) @ {speed} MB/s", end='', flush=True)
                    else:
                        dots = '.' * (int(time.time() * 2) % 4)
                        print(f"Downloading{dots} ({downloaded/1e6} MB)", end='', flush=True)
            
            print(f"Saved: {filepath} ({downloaded/1e6} MB)", flush=True)
            
    except URLError as e:
        print(f"[ERROR] Download failed: {e}", flush=True)
        # Fallback: создаём пустой файл для демо-режима
        if DEMO_MODE:
            print("[WARN] Using fallback mode for demo", flush=True)
            return False
        sys.exit(1)
    return True

def create_demo_adata(n_cells: int = 2000, n_genes: int = 500):
    print_step(f"Creating demo AnnData: {n_cells} cells × {n_genes} genes", demo=True)
    
    from scipy.sparse import random
    from scipy.stats import nbinom
    
    mean_counts = 5
    dispersion = 0.5
    
    data = random(
        n_cells, n_genes, 
        density=0.1,  
        format='csr',
        data_rvs=lambda s: nbinom.rvs(n_genes, mean_counts / (mean_counts + dispersion), size=s)
    ).astype(np.float32)
    
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['Neuron', 'Astrocyte', 'Microglia', 'Oligo'], n_cells),
        'donor_id': np.random.choice(['D1', 'D2', 'D3'], n_cells),
        'total_UMIs': np.random.poisson(5000, n_cells)
    }, index=[f"cell_{i}" for i in range(n_cells)])
    
    var = pd.DataFrame({
        'gene_name': [f"GENE_{i:04d}" for i in range(n_genes)],
        'feature_type': 'gene'
    }, index=[f"ENSG{i:011d}" for i in range(n_genes)])
    
    adata = sc.AnnData(X=data, obs=obs, var=var)
    adata.obs['sex'] = np.random.choice(['M', 'F'], n_cells)
    
    print_step(f"✓ Demo data created", demo=True)
    return adata

def setup_directories():
    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def load_or_create_data():
    if DEMO_MODE:
        print_step("Running in DEMO mode — using synthetic data", demo=True)
        return create_demo_adata(DEMO_CELLS)
    
    data_path = Path(DATA_DIR) / DATA_FILE
    
    if not data_path.exists():
        success = download_with_progress(DATA_URL, str(data_path))
        if not success:
            print("[WARN] Falling back to demo mode", flush=True)
            return create_demo_adata(DEMO_CELLS)
    
    print(f"[INFO] Loading AnnData from {data_path}...")
    adata = sc.read_h5ad(str(data_path))
    print(f"[OK] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    
    if os.getenv('QUICK_TEST', 'false').lower() == 'true':
        print("[INFO] QUICK_TEST mode: subsampling 5000 cells")
        adata = adata[:5000].copy()
    
    return adata

def run_qc(adata):
    print_step("Running QC metrics...")
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata = adata[:, adata.var["total_counts"] < 1e7]
    return adata

def remove_doublets(adata):
    if DEMO_MODE:
        print_step("Skipping Scrublet in demo mode (heavy operation)", demo=True)
        return adata
    
    print_step("Running Scrublet for doublet detection...")
    try:
        if adata.n_obs > 10000 and not DEMO_MODE:
            print("[WARN] Large dataset: consider using DEMO_MODE for faster testing")
        sce.pp.scrublet(adata)
        adata = adata[~adata.obs['predicted_doublet'], :]
        print(f"[OK] After doublet removal: {adata.n_obs} cells")
    except Exception as e:
        print(f"[WARN] Scrublet skipped: {e}")
    return adata

def filter_low_counts(adata, min_gene, min_cell):
    if DEMO_MODE:
        min_gene = max(100, min_gene // 10)
        min_cell = max(500, min_cell // 10)
        print_step(f"Using relaxed filters (demo): genes≥{min_gene}, cells≥{min_cell}", demo=True)
    
    print(f"[INFO] Filtering: genes≥{min_gene}, cells≥{min_cell}")
    
    cell_col = 'total_UMIs' if 'total_UMIs' in adata.obs else 'n_counts'
    cell_counts = adata.obs[cell_col]
    gene_counts = adata.var['total_counts']
    
    adata = adata[:, gene_counts >= min_gene]
    adata = adata[cell_counts >= min_cell, :]
    
    print(f"[OK] After filtering: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata

def run_dimensionality_reduction(adata):
    print_step("Running PCA...")
    n_pcs = 10 if DEMO_MODE else 50
    sc.tl.pca(adata, n_comps=min(n_pcs, min(adata.shape)-1))
    
    print_step("Running UMAP...")
    sc.pp.neighbors(adata, n_neighbors=10 if DEMO_MODE else 15)
    sc.tl.umap(adata)
    return adata

def save_outputs(adata, output_dir):
    print(f"[INFO] Saving results to {output_dir}...")
    
    sc.pp.log1p(adata)
    
    processed_path = Path(output_dir) / "processed_adata.h5ad"
    adata.write(processed_path)
    print(f"[OK] Saved: {processed_path}")
    
    color_cols = [c for c in ['supercluster_term', 'donor_id', 'cell_type'] if c in adata.obs.columns]
    if not color_cols:
        color_cols = ['sex'] if 'sex' in adata.obs else []
    
    for col in color_cols[:2]:
        try:
            sc.pl.umap(adata, color=col, save=f"_{col}.png", show=False)
            print(f"[OK] Plot: umap_{col}.png")
        except Exception as e:
            print(f"[WARN] Could not plot {col}: {e}")
    
    try:
        sc.pl.pca_variance_ratio(adata, n_pcs=20, save="_pca.png", show=False)
        print("[OK] Plot: pca_variance.png")
    except:
        pass
    
    report_path = Path(output_dir) / "summary.txt"
    with open(report_path, 'w') as f:
        f.write(f"scRNA-seq Analysis Summary\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Genes: {adata.n_vars}\n")
        f.write(f"Demo mode: {DEMO_MODE}\n")
        if 'cell_type' in adata.obs:
            f.write(f"\nCell types:\n{adata.obs['cell_type'].value_counts().head(5)}\n")
    print(f"[OK] Report: {report_path}")

def run_web_server(output_dir, port):
    from flask import Flask, jsonify, send_from_directory
    app = Flask(__name__)
    
    @app.route('/')
    def index():
        return jsonify({
            "service": "scRNA-seq Analysis Results",
            "status": "ready",
            "demo_mode": DEMO_MODE,
            "output_dir": output_dir
        })
    
    @app.route('/files')
    def list_files():
        files = [f.name for f in Path(output_dir).glob('*') if f.is_file()]
        return jsonify({"files": files})
    
    @app.route('/plot/<path:filename>')
    def serve_plot(filename):
        return send_from_directory(output_dir, filename)
    
    print(f"Web server: http://0.0.0.0:{port}", flush=True)
    app.run(host='0.0.0.0', port=port)

def main():
    mode_str = "DEMO" if DEMO_MODE else "Full mode"
    print(f"{mode_str} scRNA-seq Analysis Pipeline")
    print(f"Config: DATA_DIR={DATA_DIR}, OUTPUT_DIR={OUTPUT_DIR}, WEB_MODE={WEB_MODE}")
    if DEMO_MODE:
        print(f"Demo params: {DEMO_CELLS} cells, relaxed filters")
    print("-" * 50, flush=True)
    
    setup_directories()
    adata = load_or_create_data()
    
    adata = run_qc(adata)
    adata = remove_doublets(adata)
    adata = filter_low_counts(adata, MIN_COUNTS_PER_GENE, MIN_COUNTS_PER_CELL)
    adata = run_dimensionality_reduction(adata)
    save_outputs(adata, OUTPUT_DIR)
    
    print(f"{mode_str} analysis complete", flush=True)
    
    if WEB_MODE:
        run_web_server(OUTPUT_DIR, APP_PORT)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
