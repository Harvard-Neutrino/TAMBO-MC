FROM ubuntu:22.04 AS base

# Install a bunch of tools
RUN apt-get update && \
    apt-get install -yq tzdata && \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

RUN apt-get -y --no-install-recommends install \
    vim \
    build-essential \
    gcc \
    unzip \
    g++ \
    git \
    libssl-dev \
    ca-certificates \
    wget \
    libgsl-dev \
    pkg-config \
    libhdf5-serial-dev \
    libboost-all-dev \
    python3 \
    python3-pip \
    cmake \ 
    gfortran \
    git \
    rsync \
    libpthread-stubs0-dev \
    graphviz  && \
    rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2.tar.gz --no-check-certificate && \
    tar -zxvf cmake-3.22.2.tar.gz && \
    cd cmake-3.22.2 && ./bootstrap && make -j4 && make install

# Install Julia
WORKDIR /opt/
RUN wget https://github.com/JuliaLang/julia/releases/download/v1.10.5/julia-1.10.5.tar.gz --no-check-certificate && \
    tar -zxvf julia-1.10.5.tar.gz
WORKDIR /opt/julia-1.10.5
RUN make
ENV PATH="${PATH}:/opt/julia-1.10.5/"
ENV PYTHONPATH=/home/awqaykamayuq/.local/lib/

# Install conda
ENV PATH="/opt/miniconda3/bin:${PATH}"
ARG PATH="/opt/miniconda3/bin:${PATH}"
ENV PYTHONPATH=""
RUN cd /opt &&\
    wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --no-check-certificate\
    && mkdir /opt/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3\
    && rm -f Miniconda3-latest-Linux-x86_64.sh
ENV LD_LIBRARY_PATH=/opt/miniconda3/lib

# Install python
RUN conda install -y python=3.11 
RUN conda install pip
RUN python -m pip install --upgrade pip

RUN useradd -r awqaykamayuq
RUN mkdir /home/awqaykamayuq
RUN mkdir -p /home/awqaykamayuq/.local/lib

# Install TauRunner somewhere the user can modify cuz of tables
WORKDIR /
RUN wget https://github.com/icecube/TauRunner/archive/refs/heads/master.zip
RUN unzip master.zip 
WORKDIR /TauRunner-master
RUN python3 -m pip install --target /home/awqaykamayuq/.local/lib .
RUN python3 -m pip install --target /home/awqaykamayuq/.local/lib proposal

# Install FLUKA
RUN mkdir -p /home/awqaykamayuq/.local/src/fluka/
COPY ./resources/fluka2024.1-linux-gfor64bit-9.4-glibc2.17-AA.tar.gz /home/awqaykamayuq/.local/src/fluka
WORKDIR /home/awqaykamayuq/.local/src/fluka
RUN tar -xvf fluka2024.1-linux-gfor64bit-9.4-glibc2.17-AA.tar.gz 
RUN rm fluka2024.1-linux-gfor64bit-9.4-glibc2.17-AA.tar.gz 
ENV FLUPRO=/home/awqaykamayuq/.local/src/fluka/
ENV FLUFOR=gfortran
RUN make

# Install CORSIKA
RUN pip install conan
ENV CORSIKA_BUILD=/home/awqaykamayuq/.local/build/corsika/
ENV CORSIKA_INSTALL=/home/awqaykamayuq/.local/lib/
ENV CORSIKA_SRC=/home/awqaykamayuq/.local/src/corsika/
RUN mkdir -p ${CORSIKA_BUILD} ${CORSIKA_INSTALL}
RUN git clone --recursive https://gitlab.iap.kit.edu/AirShowerPhysics/corsika.git ${CORSIKA_SRC}
WORKDIR ${CORSIKA_BUILD}
RUN conda install -y -c conda-forge gcc=12.1.0
RUN conda install -y -c conda-forge ncurses
RUN $CORSIKA_SRC/conan-install.sh --source-directory $CORSIKA_SRC/ --release
RUN $CORSIKA_SRC/corsika-cmake.sh -c "-DCMAKE_BUILD_TYPE="Release" -DWITH_FLUKA=ON -DCMAKE_INSTALL_PREFIX=$CORSIKA_INSTALL"

WORKDIR /home/awqaykamayuq
RUN mkdir TAMBO-MC/ TAMBO-MC/Tambo/ TAMBO-MC/resources/

COPY ./Tambo /home/awqaykamayuq/TAMBO-MC/Tambo/
COPY ./resources /home/awqaykamayuq/TAMBO-MC/resources
COPY ./resources/proposal_tables /home/awqaykamayuq/TAMBO-MC/resources/proposal_tables
COPY ./resources/proposal_tables /home/awqaykamayuq/.local/lib/taurunner/resources/proposal_tables/
COPY ./.TAMBO.txt /home/awqaykamayuq/
ENV PYTHONPATH=/home/awqaykamayuq/.local/lib/
RUN echo "cat ~/.TAMBO.txt" > /home/awqaykamayuq/.bashrc
RUN julia -e 'using Pkg; Pkg.activate("./TAMBO-MC/Tambo"); Pkg.instantiate()'

RUN chown -R awqaykamayuq /home/awqaykamayuq
USER awqaykamayuq
SHELL ["/bin/bash", "-c"]
ENTRYPOINT [ "/bin/bash" ]
