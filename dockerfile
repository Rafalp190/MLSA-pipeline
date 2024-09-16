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
    muscle \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

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
COPY Snakefile /app/Snakefile

RUN pip3 install --no-cache-dir -r /app/requirements.txt
RUN pip3 install --no-cache-dir snakemake
# Copy the src directory into the container (contains the Python scripts)
COPY src/ /app/src/

# Set the working directory to /app
WORKDIR /app

# Specify the entry point for Snakemake
CMD ["sh", "-c", "snakemake --snakefile Snakefile --cores ${CORES:-1}"]
