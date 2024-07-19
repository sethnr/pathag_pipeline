FROM continuumio/miniconda3

# chage to root dir, add conda to path
WORKDIR /root/
ENV PATH /opt/conda/envs/env/bin:$PATH

#install conda packages from environment (nb: linux64 arch)
COPY environment.yml .
#RUN conda env update -n root -f environment.yml && conda clean -a

RUN conda install -c conda-forge numpy
RUN conda install -c conda-forge pandas
RUN conda install -c conda-forge scikit-learn
RUN conda install -c conda-forge scikit-allel
RUN conda install -c conda-forge gsl


ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"
