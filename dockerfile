# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    curl \
    build-essential \
    make \
    wget \
    clustalw \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Ensure Python is in the PATH
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install MUSCLE
RUN apt-get update && apt-get install -y wget && \
    wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64 -O /usr/bin/muscle && \
    chmod +x /usr/bin/muscle

# Install SeqKit (for concatenation)
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_amd64.tar.gz \
    && tar -zxvf seqkit_linux_amd64.tar.gz \
    && chmod +x seqkit \
    && mv seqkit /usr/local/bin/seqkit \
    && rm seqkit_linux_amd64.tar.gz

# Install MrBayes (for Bayesian phylogenetic inference)
RUN wget https://github.com/NBISweden/MrBayes/archive/v3.2.7a.tar.gz \
    && tar -xvzf v3.2.7a.tar.gz \
    && cd MrBayes-3.2.7a \
    && ./configure --enable-mpi \
    && make \
    && make install \
    && cd .. \
    && rm -rf MrBayes-3.2.7a v3.2.7a.tar.gz

# Install Python packages: BioPhylo, DendroPy, and AMAS
COPY requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir -r /app/requirements.txt
RUN pip3 install --no-cache-dir snakemake
RUN apt-get update && apt-get install -y wget && \
    wget https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 -O /usr/bin/yq && \
    chmod +x /usr/bin/yq

# Copy the Snakefile
COPY Snakefile /app/Snakefile

# Copy the src directory into the container (contains the Python scripts)
COPY src/ /app/src/

# Set the working directory to /app
WORKDIR /app

# Specify the entry point for Snakemake or a Python command for testing
#CMD ["python3", , "-h"]
CMD ["sh", "-c", "cores=$(yq e '.cores // 1' /mnt/config.yaml) && snakemake --snakefile Snakefile --configfile /mnt/config.yaml --cores $cores | tee /mnt/output/log/snakemake.log"]