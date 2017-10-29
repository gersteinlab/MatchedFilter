FROM ubuntu:14.04

RUN apt-get update && apt-get install -y \
	curl \
	bzip2 \
	git

RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /opt/miniconda.sh

RUN bash /opt/miniconda.sh -bfp /opt

RUN rm -rf /opt/miniconda.sh

RUN /opt/bin/conda install -y python=2

RUN /opt/bin/conda update conda

RUN /opt/bin/conda clean --all --yes

RUN /opt/bin/conda install -y -c bioconda metaseq

RUN /opt/bin/conda install -y seaborn

WORKDIR /opt

RUN git clone https://github.com/gersteinlab/MatchedFilter.git

COPY Dockerfile /opt/

MAINTAINER Anurag Sethi, Seven Bridges, <anurag.sethi@sbgenomics.com>
