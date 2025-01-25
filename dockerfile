# Base image with Miniconda pre-installed
FROM continuumio/miniconda3:latest AS builder

# Copy env yml file
COPY environment_cluster.yml .

# Install dependencies using mamba for faster resolution and clean up
RUN conda install -n base -c conda-forge mamba && \
    mamba env create -f environment_cluster.yml && \
    conda clean -a && \
    echo "source activate clustering_analysis_env" > ~/.bashrc

# Copy the Python module scripts
COPY py_scripts/cluster_module_files ./cluster_module_files
COPY py_scripts/main_cluster_consensus.py .

# Precompile Python scripts to improve runtime performance
RUN conda run -n clustering_analysis_env python3 -m compileall -q .

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH" \
    PYTHONWARNINGS="ignore::DeprecationWarning,ignore::FutureWarning"

# Copy the environment and precompiled scripts from the builder
WORKDIR /app
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /app /app

# Allow mounting external directories
RUN mkdir -p /data /test && chmod -R 777 /data /test

# Use ENTRYPOINT for modular CLI handling
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "clustering_analysis_env", "python3", "/app/main_cluster_consensus.py"]
