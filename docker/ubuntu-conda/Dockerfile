FROM ubuntu:18.04
MAINTAINER IDAES Tech Team <idaes-support@idaes.org>

ARG IUSER=idaes

RUN echo "\n____ INSTALL PACKAGES ___\n" \
    && apt-get -qq update \
    && apt-get -qq -y upgrade \
    && apt-get -qq -y install curl bzip2 locales \
    && apt-get -qq -y install build-essential libgfortran4 liblapack-dev \
    && apt-get -qq -y install  openssh-client git \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN echo "\n___ SET LOCALE AND TIMEZONE ___\n" \
    && echo "Etc/UTC" > /etc/timezone \
    && update-locale LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN useradd --no-log-init --create-home --shell /bin/bash $IUSER

USER $IUSER
WORKDIR /home/$IUSER
ENV PATH=/home/$IUSER/.idaes/bin:/home/$IUSER/miniconda3/bin:$PATH:/home/$IUSER/.local/bin

RUN echo "\n___ INSTALL CONDA ___\n" \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
       -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bf \
    && rm -rf /tmp/miniconda.sh

RUN echo "\n___ CREATE BASE CONDA ENV ___\n" \
    && conda install -y python=3 \
    && conda update conda \
    && conda clean --all --yes
