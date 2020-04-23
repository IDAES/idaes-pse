FROM ubuntu:18.04
MAINTAINER IDAES Tech Team <idaes-support@idaes.org>

ARG IUSER=idaes

RUN apt-get -qq update \
    && apt-get -qq -y upgrade \
    && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
       -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y install build-essential libgfortran3 openssh-client git \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && echo "Etc/UTC" > /etc/timezone \
    && conda clean --all --yes \
    && useradd --no-log-init --create-home --shell /bin/bash $IUSER

USER $IUSER
ENV PATH /opt/conda/bin:$PATH
WORKDIR /home/$IUSER
ENV PATH=$PATH:/home/$IUSER/.idaes/bin
RUN conda init bash
