#To run use:  docker build -t mach3 .
# KS: Get glorious container from Luke which will work as a base
FROM picker24/root_v6_26_10:alma9 AS mach3_build

# Declare the build argument
ARG MACH3_VERSION
ENV MACH3_VERSION=${MACH3_VERSION:-develop}

ARG CMAKE_OPTIONS
ENV CMAKE_OPTIONS=${CMAKE_OPTIONS:-DMaCh3_PYTHON_ENABLED=ON}

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
RUN cmake ${CMAKE_OPTIONS} ${MACH3_WORK_DIR}
RUN make VERBOSE=1 && make install

# KS: Need to set them here, otherwise container using this container will not be able to find MaCh3
ENV MaCh3_ROOT=${MACH3_INSTALL_DIR}
ENV PATH=${MaCh3_ROOT}/bin:${PATH} \
    LD_LIBRARY_PATH=${MaCh3_ROOT}/lib:${LD_LIBRARY_PATH}

# Start from MaCh3 install dir
WORKDIR ${MACH3_INSTALL_DIR}
