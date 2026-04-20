FROM python:3.11-slim

ENV DATA_URL=https://datasets.cellxgene.cziscience.com/b9cbe943-ad26-4cac-8798-6453b80834bf.h5ad \
    DATA_FILE=HumanBrain_NuclAccumb.h5ad \
    DATA_DIR=/data \
    OUTPUT_DIR=/output \
    MIN_COUNTS_PER_GENE=1000 \
    MIN_COUNTS_PER_CELL=15000 \
    APP_PORT=8080 \
    WEB_MODE=false \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libhdf5-dev \
    libopenblas-dev \
    liblapack-dev \
    libgomp1 \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 appuser

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY scrna_seq_data_analysis.py .

RUN mkdir -p ${DATA_DIR} ${OUTPUT_DIR} && \
    chown -R appuser:appuser /app ${DATA_DIR} ${OUTPUT_DIR}

USER appuser

EXPOSE ${APP_PORT}

CMD ["python", "scrna_seq_data_analysis.py"]
