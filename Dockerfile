# This Dockerfile is adapted (with modifications) from this Dockerfile:
# https://github.com/jupyter/docker-stacks/blob/master/scipy-notebook/Dockerfile

ARG BASE_CONTAINER=jupyter/minimal-notebook
FROM $BASE_CONTAINER

MAINTAINER Project IDAES <hhelgammal@lbl.gov>

# Maintainer Note: We're using bokeh for plotting (not matplotlib). Uncomment if matplotlib animation needed.
# USER root
# ffmpeg for matplotlib anim
#RUN apt-get update && \
#    apt-get install -y --no-install-recommends ffmpeg && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*

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
    jupyter labextension install @jupyter-widgets/jupyterlab-manager@^0.37.0 && \
    # Pin to 0.6.2 until we can move to Lab 0.35 (jupyterlab_bokeh didn't bump to 0.7.0)
    jupyter labextension install jupyterlab_bokeh@0.6.2 && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    rm -rf /home/$NB_USER/.node-gyp && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# Maintainer Note: This is an ML plotting library. Uncomment if you need it.
# Install facets which does not have a pip or conda package at the moment
# RUN cd /tmp && \
#    git clone https://github.com/PAIR-code/facets.git && \
#    cd facets && \
#    jupyter nbextension install facets-dist/ --sys-prefix && \
#    cd && \
#    rm -rf /tmp/facets && \
#    fix-permissions $CONDA_DIR && \
#    fix-permissions /home/$NB_USER

# Maintainer Note: We're using bokeh for plotting (not matplotlib). Uncomment if matplotlib needed.
# Import matplotlib the first time to build the font cache.
# ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
# RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot" && \
#    fix-permissions /home/$NB_USER

# Add idaes directory and change permissions to the notebook user:
ADD . /home/$NB_USER/idaes
USER root
RUN sudo apt-get update
# ENV DEBIAN_FRONTEND noninteractive
# RUN sudo DEBIAN_FRONTEND=noninteractive sudo apt-get install -y tzdata
RUN echo "America/Los_Angeles" > /etc/timezone
RUN chown -R $NB_UID /home/$NB_USER/idaes

# Install ipopt:
USER $NB_UID
RUN conda install -c conda-forge ipopt 

# Install idaes requirements.txt
WORKDIR /home/$NB_USER/idaes
RUN pip install -r requirements.txt

RUN python setup.py install

WORKDIR /home/$NB_USER/
USER $NB_UID