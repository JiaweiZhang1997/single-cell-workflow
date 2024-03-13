# 使用官方Ubuntu作为基础镜像
FROM ubuntu:latest

# 更新系统软件包并安装必要依赖
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        wget \
        unzip \
        libglib2.0-0 \
        zlib1g \
        libglu1-mesa-dev \
        libsm6 \
        libxrender1 \
        libxtst6 \
        default-jdk \
        python3 \
        python3-pip \
        libssl-dev \
        git \
        openssh-client && \
    rm -rf /var/lib/apt/lists/*

# install miniconda
ENV MINICONDA_VER=4.9.2
RUN curl -sSLO https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh && \
    bash Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh

# add conda to PATH
ENV PATH=/usr/local/miniconda/bin:$PATH
RUN conda update -n base -c defaults conda

# create work dir
WORKDIR /opt/cellranger

# download and install Cell Ranger
ARG CELLRANGER_VERSION=6.1.2
RUN wget https:// downloads.10xgenomics.com/supp ort/single-cell-gene-expression/software/cellranger-${CELLRANGER_VERSION}.tar.gz && \
    tar -xzvf cellranger-${CELLRANGER_VERSION}.tar.gz && \
    rm cellranger-${CELLRANGER_VERSION}.tar.gz

# 设置环境变量
ENV PATH="/opt/cellranger/cellranger-${CELLRANGER_VERSION}:${PATH}"

# 设置容器启动时执行命令
CMD ["cellranger", "--version"]

# 使用COPY指令将额外的文件或配置要复制到镜像中
# COPY my_config_file.txt /path/in/container/my_config_file.txt