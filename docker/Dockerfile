FROM nfcore/base
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing base requirements for pathway analysis"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/dean/bin:$PATH