FROM ubuntu:22.04

# ARGUMENTS
###########
ARG OPENSTRUCTURE_VERSION="2.4.0"
ARG SRC_FOLDER="/usr/local/src"
ARG CPUS_FOR_MAKE=2
ARG OPENMM_VERSION="7.7.0"
ARG OPENMM_INCLUDE_PATH="/usr/local/openmm/include/"
ARG OPENMM_LIB_PATH="/usr/local/openmm/lib/"
ARG DEBIAN_FRONTEND=noninteractive

# INSTALL SYSTEM DEPS
#####################
RUN apt-get update -y && apt-get install -y cmake \
                                            g++ \
                                            gfortran \
                                            wget \
                                            tar \
                                            libsqlite3-dev \
                                            sip-dev \
                                            libtiff-dev \
                                            libfftw3-dev \
                                            libeigen3-dev \
                                            libboost-all-dev \
                                            libpng-dev \
                                            python3-all \
                                            python3-numpy \
                                            python3-scipy \
                                            python3-pandas \
                                            python3-networkx \
                                            doxygen \
                                            swig \
                                            clustalw \
                                            cython3 \
                                            voronota \
                                            locales && \
                                            # CLEANUP
                                            rm -rf /var/lib/apt/lists/*

# INSTALL OPENMM
################
RUN cd ${SRC_FOLDER} && \
    wget -O openmm-${OPENMM_VERSION}.tar.gz -nc https://github.com/pandegroup/openmm/archive/${OPENMM_VERSION}.tar.gz && \
    mkdir ${SRC_FOLDER}/openmm-${OPENMM_VERSION} && \
    tar xf openmm-${OPENMM_VERSION}.tar.gz -C ${SRC_FOLDER}/openmm-${OPENMM_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build && \
    cd ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build && \
    cmake .. && make -j $CPUS_FOR_MAKE && make install && \
    cd ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build/python && \
    python3 setup.py build && python3 setup.py install && \
    rm ${SRC_FOLDER}/openmm-${OPENMM_VERSION}.tar.gz && \
    rm -rf ${SRC_FOLDER}/openmm-${OPENMM_VERSION}

# INSTALL OST
#############
RUN cd ${SRC_FOLDER} && \
    # copy ost release
    wget -O openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -nc https://git.scicore.unibas.ch/schwede/openstructure/-/archive/${OPENSTRUCTURE_VERSION}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz && \
    mkdir openstructure-${OPENSTRUCTURE_VERSION} && \
    tar xf openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -C ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cd ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cmake .. -DOPTIMIZE=ON \
             -DENABLE_MM=ON \
             -DCOMPILE_TMTOOLS=1 \
             -DOPEN_MM_LIBRARY=$OPENMM_LIB_PATH/libOpenMM.so \
             -DOPEN_MM_INCLUDE_DIR=$OPENMM_INCLUDE_PATH \
             -DOPEN_MM_PLUGIN_DIR=$OPENMM_LIB_PATH/plugins \
             -DENABLE_GFX=ON \
             -DENABLE_GUI=OFF \
             -DENABLE_INFO=OFF \
             -DCMAKE_C_FLAGS="-isystem /usr/include/boost/ -isystem ${OPENMM_INCLUDE_PATH}/include" \
             -DCMAKE_CXX_FLAGS="-isystem /usr/include/boost/ -isystem ${OPENMM_INCLUDE_PATH}/include" && \
    make -j ${CPUS_FOR_MAKE} && \
    wget ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz && \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb && stage/bin/chemdict_tool update ../modules/conop/data/charmm.cif compounds.chemlib charmm && \
    cmake .. -DCOMPOUND_LIB=${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build/compounds.chemlib && \
             make -j ${CPUS_FOR_MAKE} && make check && make install && \
    rm ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz && \
    rm -rf ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}

# ENVIRONMENT
#############
WORKDIR /home
ENV OST_ROOT="/usr/local"
ENV PYTHONPATH="/usr/local/lib64/python3.10/site-packages"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib64:/usr/local/openmm/lib/"
ENV OPENSTRUCTURE_VERSION=$OPENSTRUCTURE_VERSION

ENTRYPOINT ["ost"]
CMD ["-i"]

