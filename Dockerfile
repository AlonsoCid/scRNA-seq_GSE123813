# Miniconda base image
FROM continuumio/miniconda3

# Add environment file to container
COPY environment.yml /tmp/environment.yml

# Create the environment
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# Set environment PATH
ENV PATH="/opt/conda/envs/scrna/bin:$PATH"

# Copy project files
COPY scripts/ /scripts/
COPY Snakefile /Snakefile

# Set working directory
WORKDIR /data
