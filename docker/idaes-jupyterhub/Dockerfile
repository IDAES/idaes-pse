# This Dockerfile is adapted (with modifications) from this Dockerfile:
# https://github.com/jupyter/docker-stacks/blob/master/scipy-notebook/Dockerfile

ARG BASE_CONTAINER=jupyter/minimal-notebook
FROM $BASE_CONTAINER

MAINTAINER Project IDAES <hhelgammal@lbl.gov>

USER $NB_UID

# Install Python 3 packages
RUN conda install --quiet --yes 'ipywidgets' 'xlrd'  && \
    # Activate ipywidgets extension in the environment that runs the notebook server
    jupyter nbextension enable --py widgetsnbextension --sys-prefix && \
    # Also activate ipywidgets extension for JupyterLab
    # Check this URL for most recent compatibilities
    # https://github.com/jupyter-widgets/ipywidgets/tree/master/packages/jupyterlab-manager
    jupyter labextension install @jupyter-widgets/jupyterlab-manager && \
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

# Add top-level idaes source directory and change permissions to the notebook user:
ADD . /home/idaes
USER root
RUN sudo apt-get update
RUN sudo apt-get -y install libgfortran3
RUN echo "America/Los_Angeles" > /etc/timezone
RUN chown -R $NB_UID /home/idaes

# Install idaes requirements.txt
USER $NB_UID
WORKDIR /home/idaes
RUN pip install -r requirements-dev.txt
RUN python setup.py install
RUN idaes get-extensions

WORKDIR /home
ENV PATH=$PATH:/home/jovyan/.idaes/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jovyan/.idaes/lib:/opt/conda/lib/
USER $NB_UID