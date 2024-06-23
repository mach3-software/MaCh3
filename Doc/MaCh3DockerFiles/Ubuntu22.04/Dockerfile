#To run use:  docker build -t mach3 .
# KS: Get glorious container from ROOT which will work as a base
FROM rootproject/root:6.32.02-ubuntu22.04 AS mach3_build

RUN apt update && apt upgrade -y

RUN apt-get install -y nlohmann-json3-dev
#MISC_SW
RUN apt install -y --no-install-recommends vim less nano gdb csh tcsh ed quota python3 python3-dev cvs procmail ca-certificates

ENV MACH3_WORK_DIR=/opt/MaCh3/
ENV MACH3_INSTALL_DIR=/opt/MaCh3/build/

RUN mkdir -p ${MACH3_WORK_DIR}

WORKDIR /opt/
# KS: Let's clone MaCh3
RUN --mount=type=ssh git clone https://github.com/mach3-software/MaCh3 ${MACH3_WORK_DIR}
WORKDIR ${MACH3_WORK_DIR}
RUN git checkout ${MACH3_VERSION}

RUN mkdir -p ${MACH3_INSTALL_DIR}
WORKDIR ${MACH3_INSTALL_DIR}
RUN cmake ${MACH3_WORK_DIR}
RUN make VERBOSE=1
RUN make install

# KS: Need to set them here, otherwise container using this container will not be able to find MaCh3
ENV MaCh3_ROOT=${MACH3_INSTALL_DIR}
ENV PATH=${MaCh3_ROOT}/bin:${PATH} \
    LD_LIBRARY_PATH=${MaCh3_ROOT}/lib:${LD_LIBRARY_PATH}

# Start from MaCh3 install dir
WORKDIR ${MACH3_INSTALL_DIR}
