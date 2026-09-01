FROM python:3.12-slim

WORKDIR /app

# Install system dependencies for lrge
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install lrge
RUN curl -L https://github.com/mbhall88/lrge/releases/latest/download/lrge-x86_64-unknown-linux-musl.tar.gz \
    | tar -xz -C /usr/local/bin && chmod +x /usr/local/bin/lrge

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy script
COPY manifest_generator.py .

ENTRYPOINT ["python", "manifest_generator.py"]
