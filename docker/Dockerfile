FROM vinuesa/get_phylomarkers:latest

LABEL authors="Pablo Vinuesa, Bruno Contreras Moreira"
LABEL description="Container with GET_HOMOLOGUES⁠ and GET_PHYLOMARKERS⁠ ready  use."
LABEL release.notes="GET_HOMOLOGUES ships with Pfam-A and swissprot installed."
LABEL license="GPLv3"
LABEL maintainers="bcontreras@eead.csic.es, vinuesa@ccg.unam.mx"

# install unsatisfied get_homologues dependencies
USER root
RUN apt update && apt install --no-install-recommends -y \ 
	libgd-gd2-perl \
	r-cran-cluster \
	r-cran-dendextend \
	r-cran-factoextra \
	&& apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN cpanm Inline::C Inline::CPP

# clone get_homologues & get executable binaries & db
WORKDIR /
RUN cachebuster=b953b36 git clone https://github.com/eead-csic-compbio/get_homologues.git
WORKDIR /get_homologues
RUN git pull
RUN chmod a+wr db
RUN make install_docker
RUN make testR 
RUN make clean cleanR
RUN chgrp -R root bin && chmod -R a+rx bin

## <---- remove some uneeded folders to save space 
RUN rm -rf .git sample_intersection sample_*homologues
RUN rm -rf /get_plylomarkers/.git
WORKDIR /get_plylomarkers/bin
RUN rm -rf .git clustal-omega-1.2.4 FastTree_v2.1.11 macosx-intel parallel-20171022 PhiPack phylip-3.695

# add version name to image
ARG version
LABEL version=$version
RUN echo $version

# prepare env
WORKDIR /home/you
ENV PATH="/get_homologues:${PATH}"
USER you
