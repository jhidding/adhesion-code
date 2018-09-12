FROM debian:testing

# == Dependencies ==
RUN apt-get update && apt-get install -y --no-install-recommends \
    argagg-dev \
    cmake \
    g++ \
    git \
    haskell-stack \
    libcgal-dev \
    libfftw3-dev \
    libfmt-dev \
    libhdf5-dev \
    libyaml-cpp-dev \
    lmodern \
    make \
    pkg-config \
    rsync \
    texlive \
    texlive-fonts-extra \
    texlive-latex-extra \
    texlive-latex-recommended \
    texlive-xetex \
    wget

ENV PATH /root/.local/bin:$PATH

# == We need at least version 2.2.3 of Pandoc ==
RUN wget https://github.com/jgm/pandoc/releases/download/2.2.3.2/pandoc-2.2.3.2-1-amd64.deb
RUN dpkg -i ./pandoc-2.2.3.2-1-amd64.deb

# == We need the latest version of XTensor ==
WORKDIR /root
RUN git clone https://github.com/QuantStack/xtl.git \
 && cd xtl \
 && mkdir build \
 && cd build \
 && cmake .. -DCMAKE_INSTALL_PREFIX=~/.local \
 && make install
RUN git clone https://github.com/QuantStack/xtensor.git \
 && cd xtensor \
 && mkdir build \
 && cd build \
 && cmake .. -DCMAKE_INSTALL_PREFIX=~/.local \
 && make install

# == Build the project and the report ==
COPY . /app
WORKDIR /app
RUN  make && make report

CMD /bin/bash
