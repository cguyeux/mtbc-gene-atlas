FROM python:3.11-slim

WORKDIR /app

COPY pyproject.toml /app/
RUN pip install --no-cache-dir -e ".[prod]"

COPY backend /app/backend
COPY frontend /app/frontend
COPY content /app/content

ENV PYTHONUNBUFFERED=1
EXPOSE 8000

# DB is populated on startup (auto_ingest) from content/gene_xref.tsv + content/genes/.
CMD ["gunicorn", "backend.app:app", "-k", "uvicorn.workers.UvicornWorker", "-b", "0.0.0.0:8000", "-w", "2"]
