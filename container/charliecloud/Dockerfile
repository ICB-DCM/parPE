FROM ubuntu:24.04

COPY /container_files /container-files

RUN /container-files/install.sh && rm -rf /tmp && mkdir /tmp

ENV BASH_ENV "/etc/drydock/.env"
ENV CC clang
ENV CXX clang++
ENV OMPI_CC clang
ENV OMPI_CXX clang++
RUN /container-files/install_parpe.sh && echo "source /parPE/build/venv/bin/activate" >> /etc/bash.bashrc

ENV PARPE_DIR "/parPE"

RUN chmod -R ugo+rwX $PARPE_DIR
