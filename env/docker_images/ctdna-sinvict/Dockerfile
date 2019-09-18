# =====================================================================
#              erikwaskiewicz/ctdna-sinvict Dockerfile
#
# For use in ctDNA pipeline project - image for running the SiNVICT
# variant caller
# <https://github.com/erikwaskiewicz/ctdna-project>
#
# =====================================================================

FROM ubuntu:14.04 
MAINTAINER Erik Waskiewicz

# Set the working directory & copy across current directory
WORKDIR /home
COPY . /home

# Install requirements
RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    git \
    software-properties-common \
    cmake \
    libncurses-dev \
    zlib1g-dev \
    patch \
    libbz2-dev \
    liblzma-dev \
# remove apt cache after running to reduce image size
 && rm -rf /var/lib/apt/lists/*

# Install sinvict and bam-readcount
RUN git clone https://github.com/sfu-compbio/sinvict.git \
 && make -C /home/sinvict \
 && git clone https://github.com/samtools/htslib.git \
 && git clone https://github.com/samtools/bcftools.git \
 && make -C /home/bcftools \
 && git clone https://github.com/genome/bam-readcount.git \
 && cd /home/bam-readcount/ \
 && cmake /home/bam-readcount/ \
 && make \
 && export PATH=$PATH:/home/bam-readcount/bin:/home/sinvict

