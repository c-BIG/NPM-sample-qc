FROM ubuntu:22.04

LABEL base_image="ubuntu:22.04"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENV PATH="/usr/local/conda/envs/npm-sample-qc/bin/:${PATH}"
ENV PATH="/usr/local/conda/envs/npm-sample-qc/share/picard-2.27.0-0/:${PATH}"
ENV PATH="/usr/local/conda/envs/npm-sample-qc/share/verifybamid2-2.0.1-7/:${PATH}"
ENV PATH="/usr/local/conda/bin/:${PATH}"

ADD conda_env.yml /conda_env.yml

RUN export DEBIAN_FRONTEND="noninteractive" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      libtool \
      libarchive-dev \
      make cmake m4 automake \
      gcc g++ \
      python3 python-pip python-dev-is-python3 \
      git \
      wget \
      gnuplot \
      datamash \
      graphviz \
      vim \
      && rm -rf /var/lib/apt/lists/*

# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/conda
RUN PATH="/usr/local/conda/bin:${PATH}"
RUN rm Miniconda3-latest-Linux-x86_64.sh
RUN conda update -n base -c defaults conda

# install AWS CLI
RUN conda install -c conda-forge -y awscli

# create conda environment from template
RUN conda env create -f /conda_env.yml
RUN conda clean -a

# install MultiQC and MultiQC_NPM plugin
# RUN git clone --depth 1 https://github.com/c-BIG/NPM-sample-qc.git

RUN git clone -b feature/vcf-qc-implementation https://github.com/c-BIG/NPM-sample-qc.git
# RUN git clone https://github.com/c-BIG/NPM-sample-qc.git
WORKDIR /NPM-sample-qc/bin/multiqc_plugins

# RUN git checkout feature/rework-aln-variant-metrics-json
RUN python3 setup.py install
# RUN python -m pip install .
# RUN pip3 install numpy scipy argparse matplotlib
