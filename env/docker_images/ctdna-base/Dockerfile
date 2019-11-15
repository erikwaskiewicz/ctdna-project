# =====================================================================
#              erikwaskiewicz/ctdna-base Dockerfile
#
# For use in ctDNA pipeline project - base image to be used by 
# pipeline if a specific image isn't being used
# <https://github.com/erikwaskiewicz/ctdna-project>
#
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda3:4.7.10
LABEL maintainer="Erik Waskiewicz"

# Set the working directory & copy across current directory
WORKDIR /var/app
COPY . /var/app
ENV PYTHONDONTWRITEBYTECODE=true

# Install awk and ps for use with Nextflow reports
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
        gawk \
        procps \
# Remove apt cache after running to reduce image size
 && rm -rf /var/lib/apt/lists/*

# Make the conda environment
RUN conda install --yes --freeze-installed \
        python=3.7 \
        bioconda::gatk4=4.1.4.0 \
        bioconda::picard=2.20.6 \
        bioconda::samtools=1.9 \
        conda-forge::pandas=0.25.1 \
        conda-forge::numpy=1.17.2 \
# Clean up the conda env
 && conda clean -afy \
 && find /opt/conda/ -follow -type f -name '*.a' -delete \
 && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
 && find /opt/conda/ -follow -type f -name '*.js.map' -delete
