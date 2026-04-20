# scRNA-seq Data Analysis

## Что делает проект
Автоматизированный пайплайн для анализа данных single-cell РНК-секвенирования (scRNA-seq): загружает данные в формате H5AD, проводит QC, фильтрацию, удаление дублетов (Scrublet), нормализацию, PCA и UMAP. 
Результаты (обработанный AnnData-объект и визуализации) сохраняются в настраиваемую директорию. 
Поддерживает веб-интерфейс.

## Как запустить

### Базовый запуск
```bash
git clone https://github.com/poandrej23-netizen/scRNA-seq_data_analysis.git
cd scRNA-seq_data_analysis
```
### Сборка образа
docker build -t scrna-analysis .

### Запуск анализа (вывод в терминал)
docker run --rm scrna-analysis

### Запуск с включённым веб-интерфейсом
docker run -d \
  -p 8080:8080 \
  -e WEB_MODE=true \
  -e APP_PORT=8080 \
  -v app-data:/data \
  -v app-output:/output \
  --name scrna-web \
  scrna-analysis

#### Открыть в браузере:
##### → http://localhost:8080
##### → http://localhost:8080/files (список файлов)
##### → http://localhost:8080/plot/umap_supercluster_term.png (график)

### Просмотр логов в реальном времени
docker logs -f scrna-web

### Остановка и удаление
docker stop scrna-web && docker rm scrna-web

### Доказательство сохранности данных

#### Первый запуск — создаём данные
docker run --rm \
  -v app-data:/data \
  -v app-output:/output \
  --name run1 \
  scrna-analysis

#### Проверка данных
docker run --rm -v app-output:/output alpine ls -la /output
##### Ожидается: processed_adata.h5ad, umap_*.png, pca_*.png

##### Перезапуск 
docker run --rm \
  -v app-data:/data \
  -v app-output:/output \
  --name run2 \
  scrna-analysis

docker run --rm -v app-output:/output alpine ls -la /output

## Какой вариант задания выбран и почему
Вариант 3, потому что веб-интерфейс удобен для просмотра графиков UMAP, PCA. Сохранение данных удобно для проведения анализа с разными датасетами и их сравнения между собой и проверки воспроизводимости.




