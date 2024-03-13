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
ENV MINICONDA_VER="latest"

RUN curl -sSLO https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh && \
    bash Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh

# add conda to PATH
ENV PATH=/usr/local/miniconda/bin:$PATH
RUN conda update -n base -c defaults conda

# create work dir
WORKDIR /opt/cellranger

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# 添加特定的Anaconda服务器镜像源以加速下载
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/

# download and install Cell Ranger
ARG CELLRANGER_VERSION=7.2.0
RUN wget -O cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1710352364&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=LRLzaJH1XmSNYi3ZSZI4HoYiAbPO3K-l6JdTbv8oEIyLspg1ft2vZTGG64YAHB6RgPq9ynMUQrrkliAgyHoktB0jhJUeeWrlcmbbe-g4D4figj-Pn74-vdYWXYc9G09FOX3MCOFXYQaE8yj3PXfxON0QTGABkD0Q15W2rmnUbDf~8p9ZI5ocNsEfeWCdmTnj5bfqAsiUpfQ-N35VqO5pgHJ3JOsKrisuDdy0ZQD570zXRdLQsbJjGAFX5zJdKUvBvI6JKcCLN7j9mbWfxZ6T3j0jru0tJUgR7uqHUKhdy6aYiNSrOnlytbD9ICJ9mk4WSA965queZDGCsULOmFqoOA__" && \
    tar -xzvf cellranger-${CELLRANGER_VERSION}.tar.gz && \
    rm cellranger-${CELLRANGER_VERSION}.tar.gz

# 设置环境变量
ENV PATH="/opt/cellranger/cellranger-${CELLRANGER_VERSION}:${PATH}"

# 设置容器启动时执行命令
CMD ["cellranger", "--version"]

# 使用COPY指令将额外的文件或配置要复制到镜像中
# COPY my_config_file.txt /path/in/container/my_config_file.txt