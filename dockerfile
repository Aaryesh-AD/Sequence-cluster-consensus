# Base image with Python and Miniconda
FROM continuumio/miniconda3:latest

# Set environment variables for Conda
ENV PATH="/opt/conda/bin:$PATH"
ENV CONDA_ENV_NAME="clustering_analysis_env"

# Create a working directory
WORKDIR /app

# Copy the environment file first to leverage Docker caching
COPY environment_cluster.yml /app/environment_cluster.yml

# Install mamba for faster environment setup and create the environment
RUN conda install -n base -c conda-forge mamba && \
    mamba env create -f /app/environment_cluster.yml && \
    conda clean -afy && \
    echo "source activate ${CONDA_ENV_NAME}" > ~/.bashrc

# Use mamba in the SHELL to speed up installations
SHELL ["conda", "run", "-n", "clustering_analysis_env", "/bin/bash", "-c"]

# Install Clustal Omega (Clustalo)
RUN apt-get update && apt-get install -y clustalo && apt-get clean

# Copy the module files and main script
COPY py_scripts/cluster_module_files /app/cluster_module_files
COPY py_scripts/main_cluster_consensus.py /app/main_cluster_consensus.py

# Precompile Python scripts to improve runtime performance
RUN python3 -m compileall -q /app

# Set the entrypoint to directly run the main script
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "clustering_analysis_env", "python3", "/app/main_cluster_consensus.py"]
