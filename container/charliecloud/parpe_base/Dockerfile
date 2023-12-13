FROM ubuntu:22.04

COPY . /container-files

RUN /container-files/install.sh && rm -rf /tmp && mkdir /tmp

ENV BASH_ENV "/etc/drydock/.env"
ENV CC clang
ENV CXX clang++
ENV OMPI_CC clang
ENV OMPI_CXX clang++

RUN /container-files/install_parpe.sh

ENV PARPE_DIR "/parPE"

RUN chmod -R ugo+rwX $PARPE_DIR
