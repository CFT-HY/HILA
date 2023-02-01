FROM ubuntu:22.04
LABEL Description="HILA environment"

ENV LLVM_VERSION 15
SHELL ["/bin/bash", "-c"]

#installing base dependencies
RUN apt-get -y update && apt-get -y install build-essential \
            libopenmpi-dev \
            libfftw3-dev \
            libomp-dev \
            git \
            make

#installing clang dependencies for hilapp
# RUN apt-get -y update && apt-get -y install clang-$LLVM_VERSION \
#                        llvm-$LLVM_VERSION \
#                        clang-tools-$LLVM_VERSION \
#                        libclang-common-$LLVM_VERSION-dev \
#                        libclang-cpp$LLVM_VERSION-dev \
#                        libclang-$LLVM_VERSION-dev \
#                        clang-format-$LLVM_VERSION

#installing clang dependencies for hilapp
RUN apt-get -y update && apt-get -y install clang-$LLVM_VERSION \
                       libclang-$LLVM_VERSION-dev

#installing hilapp
RUN cd / && \
    git clone https://github.com/CFT-HY/HILA && \
    cd /HILA/hilapp && \
    make -j4 && make install

