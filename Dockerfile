FROM debian:testing

# == Dependencies ==
RUN apt-get update && apt-get install -y --no-install-recommends \
    argagg-dev \
    cmake \
    g++ \
    git \
    libcgal-dev \
    libfftw3-dev \
    libfmt-dev \
    libhdf5-dev \
    libyaml-cpp-dev \
    lmodern \
    make \
    pkg-config \
    python3-pip \
    rsync \
    texlive \
    texlive-fonts-extra \
    texlive-latex-extra \
    texlive-latex-recommended \
    texlive-xetex \
    wget

# == Python packages ==
RUN apt-get install -y --no-install-recommends \
    python3-dev python3-wheel python3-setuptools
RUN pip3 install pandoc-eqnos pandoc-fignos

# == We need at least version 2.2.3 of Pandoc ==
RUN wget https://github.com/jgm/pandoc/releases/download/2.4/pandoc-2.4-1-amd64.deb
RUN dpkg -i ./pandoc-2.4-1-amd64.deb

RUN apt-get install -y --no-install-recommends \
    libgsl-dev

# == Build the project and the report ==
COPY . /app
WORKDIR /app
RUN  make && make html

CMD /bin/bash
