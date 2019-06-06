FROM kaitj/neurobeer:base

MAINTAINER <tkai@uwo.ca>

run mkdir -p /src/install_scripts
COPY ./install_scripts/*.sh /src/install_scripts/

ENV DEBIAN_FRONTEND noninteractive
3 RUN bash /src/install_scripts/00_install_basics.sh > log_00_basics
RUN bash /src/install_scripts/05_install_neurobeer.sh > log_05_neurobeer
# RUN bash /src/install_scripts/06_install_camino_source.sh > log_06_camino
# RUN bash /src/install_scripts/07_install_octave.sh > log_07_octave
# RUN bash /src/install_scripts/08_install_vasst_dev.sh > log_08_vasst_dev
# RUN bash /src/install_scripts/09_install_fsl.sh > log_09_fsl

# Remove install scripts
RUN rm -rf /src

# Camino
# ENV CAMINO_DIR=/opt/camino
# ENV PATH $CAMINO_DIR/bin:$PATH
# ENV LD_LIBRARY_PATH $CAMINO_DIR/lib:$LD_LIBRARY_PATH
# ENV MANPATH $CAMINO_DIR/man:$MANPATH
# ENV CAMINO_HEAP_SIZE 16000

# Vasst-dev
# ENV VASST_DEV_HOME /opt/git/vasst-dev
# ENV PIPELINE_TOOL_DIR $VASST_DEV_HOME/tools
# ENV PATH $PIPELINE_TOOL_DIR:$PATH
