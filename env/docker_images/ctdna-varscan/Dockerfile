# =====================================================================
#  erikwaskiewicz/ctdna-varscan Dockerfile
#
# For use in ctDNA pipeline project - image for running the VarScan2
# variant caller
# <https://github.com/erikwaskiewicz/ctdna-project>
#
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda3:4.7.10
LABEL maintainer="Erik Waskiewicz"

# Set the working directory & copy across current directory
# Set python so that it doesn't create .pyc files
WORKDIR /var/app
COPY . /var/app
ENV PYTHONDONTWRITEBYTECODE=true

# Install awk and ps for use with Nextflow reports
# remove apt cache after running to reduce image size
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
        gawk \
        procps \
 && rm -rf /var/lib/apt/lists/*

# Update conda, make the conda environment, clean unnecessary files afterwards
RUN conda install --yes --freeze-installed \
        python=3.7 \
        bioconda::varscan=2.4.3 \
        bioconda::samtools=1.9 \
 && conda clean -afy \
 && find /opt/conda/ -follow -type f -name '*.a' -delete \
 && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
 && find /opt/conda/ -follow -type f -name '*.js.map' -delete
