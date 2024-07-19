FROM continuumio/miniconda3

# chage to root dir, add conda to path
WORKDIR /root/
ENV PATH /opt/conda/envs/env/bin:$PATH

#install conda packages from environment (nb: linux64 arch)
COPY environment.yml .
RUN conda env update -n root -f environment.yml && conda clean -a

RUN conda install -c conda-forge python=3.9
RUN conda install -c conda-forge mash
RUN conda install -c bioconda bcftools
RUN conda install -c conda-forge gsl
RUN conda install -c bioconda cutadapt


ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"
