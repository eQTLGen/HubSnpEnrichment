FROM nfcore/base
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for hub variant enrichment analyses"

COPY environment.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/hubsnpenrichment/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
COPY temp/IGUtilityPackage_0.2.98.tar.gz /
RUN R -e "install.packages('IGUtilityPackage_0.2.98.tar.gz', repos = NULL, type = 'source', dependencies = TRUE)"
