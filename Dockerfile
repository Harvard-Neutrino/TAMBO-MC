FROM ubuntu:18.04 AS base

# Install a bunch of tools
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    build-essential \
    gcc \
    unzip \
    g++ \
    git \
    libssl-dev \
    python3-dev \
    ca-certificates \
    wget libgsl-dev pkg-config libhdf5-serial-dev libboost-all-dev python-dev && \
    # apt-get install libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2.tar.gz --no-check-certificate && \
    tar -zxvf cmake-3.22.2.tar.gz && \
    cd cmake-3.22.2 && ./bootstrap && make -j4 && make install


FROM base AS python


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
RUN conda install python=3.11
RUN conda install pip
RUN conda install poetry
RUN python -m pip install --upgrade pip

# Install TauRunner somewhere the user can modify cuz of tables
RUN useradd -r myuser          # no specific user ID
RUN mkdir /home/myuser
RUN cd /home/myuser cp /opt/.bashrc .bashrc
RUN mkdir /home/myuser/.local/ /home/myuser/.local/lib

WORKDIR /
RUN wget https://github.com/icecube/TauRunner/archive/refs/heads/master.zip
RUN unzip master.zip 
WORKDIR /TauRunner-master
RUN pip install --target /home/myuser/.local/lib .

RUN pip install --target /home/myuser/.local/lib proposal

# Download Julia
WORKDIR /opt/
RUN wget https://github.com/JuliaLang/julia/releases/download/v1.10.5/julia-1.10.5.tar.gz --no-check-certificate && \
    tar -zxvf julia-1.10.5.tar.gz
WORKDIR /opt/julia-1.10.5
RUN make
ENV PATH="${PATH}:/opt/julia-1.10.5/"
ENV PYTHONPATH=/home/myuser/.local/lib/

WORKDIR /home/myuser
RUN mkdir TAMBO-MC/ TAMBO-MC/Tambo/ TAMBO-MC/resources/

COPY ./Tambo /home/myuser/TAMBO-MC/Tambo/
COPY ./resources /home/myuser/TAMBO-MC/resources

RUN chown -R myuser /home/myuser
USER myuser
SHELL ["/bin/bash", "-c"]
ENTRYPOINT [ "/bin/bash" ]
