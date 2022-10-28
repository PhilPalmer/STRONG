FROM continuumio/miniconda3:4.9.2
LABEL authors="Phil Palmer" \
    description="Docker image containing dependencies for Strain Resolution ON Graphs https://github.com/chrisquince/STRONG"

# Install depdendencies & deep clean the apt cache to reduce image/layer size
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y libbz2-dev libreadline-dev cmake g++ zlib1g zlib1g-dev && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Install STRONG
ENV SPATH=/
RUN cd $SPATH && \
    git clone --recurse-submodules https://github.com/chrisquince/STRONG.git && \
    cd STRONG && \
    ./install_STRONG.sh

# Add STRONG and conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /STRONG/bin:/opt/conda/envs/STRONG/bin:$PATH