FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

# Set timezone
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Rome

# Update system
RUN apt-get update && apt-get upgrade -y

# Install dependencies
RUN apt-get install -y build-essential git libftdi-dev libftdi1 doxygen python3-pip libsdl2-dev curl cmake libusb-1.0-0-dev scons gtkwave libsndfile1-dev rsync python3-six \
    autoconf automake texinfo libtool pkg-config libsdl2-ttf-dev autoconf automake autotools-dev curl libmpc-dev libmpfr-dev libgmp-dev gawk build-essential bison flex texinfo gperf libtool patchutils bc zlib1g-dev \
    vim curl sudo && \
    pip3 install argcomplete pyelftools

# Install PULP-RISCV toolchain
ENV PATH="/opt/riscv/bin:${PATH}"
RUN git clone --recursive https://github.com/pulp-platform/pulp-riscv-gnu-toolchain && \
    cd pulp-riscv-gnu-toolchain && \
    ./configure --prefix=/opt/riscv --with-arch=rv32imc --with-cmodel=medlow --enable-multilib && \
    make && cd .. && rm -r pulp-riscv-gnu-toolchain 

# Install PULP-SDK
RUN mkdir /pulp
ENV PULP_RISCV_GCC_TOOLCHAIN="/opt/riscv"
RUN curl -LO https://github.com/pulp-platform/pulp-sdk/archive/refs/tags/2022.07.22.tar.gz && \
    tar -zxvf 2022.07.22.tar.gz && rm 2022.07.22.tar.gz && mv pulp-sdk-2022.07.22 /pulp/pulp-sdk && \
    cd /pulp/pulp-sdk && source ./configs/pulp-open.sh && \
    make build

# Add personalization to Vim
RUN { echo "filetype plugin indent on"; echo "set tabstop=4"; echo "set shiftwidth=4"; echo "set expandtab"; } >> .vimrc

# Create volume
RUN mkdir /fast-ica
WORKDIR /fast-ica
VOLUME /fast-ica

CMD source /pulp/pulp-sdk/configs/pulp-open.sh && exec /bin/bash
