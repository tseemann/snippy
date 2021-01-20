FROM continuumio/miniconda3

LABEL authors="Torsten Seemann" \
      description="Docker image containing dependencies for tseemann/snippy."

# Store version dependencies
ENV VER_VT="2015.11.10"
ENV VER_SNIPPY="4.6.0"

# Create conda environment
RUN conda install -c conda-forge mamba=0.6.4 \
    && mamba create --name snippy \
    && mamba install --name snippy -c bioconda/label/cf201901 vt=${VER_VT} \
    && mamba install --name snippy -c conda-forge -c bioconda -c anaconda -c defaults snippy=${VER_SNIPPY} \
    && mamba clean -a \
    && echo "source /opt/conda/bin/activate snippy" >> ~/.bashrc;

# Add conda installation directory to path
ENV PATH /opt/conda/envs/snippy/bin:$PATH
ENV IN_DOCKER_CONTAINER Yes

# Capture conda env details
RUN conda env export --name snippy > snippy.yaml
