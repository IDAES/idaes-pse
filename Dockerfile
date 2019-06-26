# This Dockerfile is adapted (with modifications) from this Dockerfile:
# https://github.com/jupyter/docker-stacks/blob/master/scipy-notebook/Dockerfile

ARG BASE_CONTAINER=jupyter/minimal-notebook
FROM $BASE_CONTAINER

MAINTAINER Project IDAES <hhelgammal@lbl.gov>

USER $NB_UID

# Install Python 3 packages
RUN conda install --quiet --yes \
    'conda-forge::blas=*=openblas' \
    'ipywidgets=7.2*' \
    'xlrd'  && \
    # Activate ipywidgets extension in the environment that runs the notebook server
    jupyter nbextension enable --py widgetsnbextension --sys-prefix && \
    # Also activate ipywidgets extension for JupyterLab
    # Check this URL for most recent compatibilities
    # https://github.com/jupyter-widgets/ipywidgets/tree/master/packages/jupyterlab-manager
    jupyter labextension install @jupyter-widgets/jupyterlab-manager@^0.38.1 && \
    # Pin to 0.6.2 until we can move to Lab 0.35 (jupyterlab_bokeh didn't bump to 0.7.0)
    jupyter labextension install jupyterlab_bokeh@0.6.3 && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    rm -rf /home/$NB_USER/.node-gyp && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# Maintainer Note: We're using bokeh for plotting (not matplotlib). Uncomment if matplotlib needed.
# Import matplotlib the first time to build the font cache.
# ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
# RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot" && \
#    fix-permissions /home/$NB_USER

# Add idaes directory and change permissions to the notebook user:
ADD . /home/idaes
USER root
RUN wget https://idaes-files.s3.amazonaws.com/public/idaes-coinbinary-static-build.zip -P /usr/local
RUN sudo apt-get update
RUN echo "America/Los_Angeles" > /etc/timezone
RUN chown -R $NB_UID /home/idaes

# Expand idaes-coinbinary contents:
WORKDIR /usr/local
RUN unzip idaes-coinbinary-static-build.zip -d /usr/local/coinor-optimization-suite-1.8
RUN cp /usr/local/coinor-optimization-suite-1.8/bin/* /usr/local/bin

# Copying part of install-solvers here:
WORKDIR /home/idaes
RUN sudo apt-get update && sudo apt-get install -y libboost-dev
RUN sudo apt-get install -y libgfortran3
RUN wget https://ampl.com/netlib/ampl/solvers.tgz
RUN tar -xf solvers.tgz
WORKDIR /home/idaes/solvers 
RUN ./configure && make
ENV ASL_BUILD=/home/idaes/solvers/sys.x86_64.Linux

# Install ipopt from conda (pinned version to avoid JupyterHub issues):
# Commenting this in favor of using COIN-OR suite binary.
# RUN conda install -c conda-forge ipopt=3.12.12=hc6e8484_1002

# Install idaes requirements.txt
USER $NB_UID
WORKDIR /home/idaes
RUN pip install -r requirements-dev.txt
RUN make
RUN python setup.py install

WORKDIR /home
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/coinor-optimization-suite-1.8/lib:/opt/conda/lib/
ENV PATH=$PATH:/usr/local/bin
# Command to smoke-test ipopt install in the Docker container:
# cd /usr/local/idaes-coinbinary-1.8.0 && ./bin/ipopt ./test/mytoy.nl -AMPL linear_solver=ma57
USER $NB_UID