# scRNA-seq Data Analysis

## What the project does
Automated pipeline for single-cell RNA sequencing (scRNA-seq) data analysis: loads data in H5AD format, performs QC, filtering, duplicate removal (Scrublet), normalization, PCA, and UMAP. 
Results (processed AnnData object and visualizations) are saved in a custom directory. 
Supports a web interface.

By default, human brain nucleus accumbens is used as real data. You can replace it with any of your data by replacing the link in line 25 of the executable .py file

## How to use

### Quick start
```bash
git clone https://github.com/poandrej23-netizen/scRNA-seq_data_analysis.git
cd scRNA-seq_data_analysis
```
### Docker build
```bash
docker build -t scrna-analysis .
```
### Analysis start
```bash
docker run --rm scrna-analysis
```
#### Demo mode (generates an artifical dataset)
```bash
docker run --rm \
  -e DEMO_MODE=true \
  -e DEMO_CELLS=2000 \
  -v app-output:/output \
  scrna-analysis
```
#### Real data subsample analysis
```bash
docker run --rm \
  -e QUICK_TEST=true \        # it uses a real file, but will take 5000 cells
  -v app-/data \
  -v app-output:/output \
  scrna-analysis
```
### Launching with the web interface enabled
```bash
docker run -d \
  -p 8080:8080 \
  -e WEB_MODE=true \
  -e APP_PORT=8080 \
  -v app-data:/data \
  -v app-output:/output \
  --name scrna-web \
  scrna-analysis
```
#### Launching with the web interface enabled in demo mode
```bash
docker run -d \
  -p 8080:8080 \
  -e DEMO_MODE=true \
  -e WEB_MODE=true \
  -v app-output:/output \
  --name scrna-demo \
  scrna-analysis
```
#### Open in a browser:
##### http://localhost:8080
##### http://localhost:8080/files (список файлов)

### Viewing logs in real time
```bash
docker logs -f scrna-web
```
### Stopping and deleting
```bash
docker stop scrna-web && docker rm scrna-web
```
### Proof of data security

#### First launch — creating data
```bash
docker run --rm \
  -v app-data:/data \
  -v app-output:/output \
  --name run1 \
  scrna-analysis
```
#### Data verification
```bash
docker run --rm -v app-output:/output alpine ls -la /output
```
##### Expected: processed_adata.h5ad, umap_*.png, pca_*.png

##### Restart
```bash
docker run --rm \
  -v app-data:/data \
  -v app-output:/output \
  --name run2 \
  scrna-analysis

docker run --rm -v app-output:/output alpine ls -la /output
```





