#To run use:  docker build -t mach3 .
# KS: Get glorious container from Kamil which will work as a base
FROM kamilskwarczynski/nukamil:latest AS mach3_build

# Add a label for the author
LABEL maintainer="The MaCh3 Collaboration"
LABEL website="https://mach3-software.github.io/MaCh3/"
LABEL compiler="GNU 11.4.1"
LABEL root_version="v6.30.04"
LABEL cuda_version="v12.6.2"
LABEL org.opencontainers.image.description="Official MaCh3 container"

# Declare the build argument
ARG MACH3_VERSION
ENV MACH3_VERSION=${MACH3_VERSION:-develop}

ARG CMAKE_OPTIONS
ENV CMAKE_OPTIONS=${CMAKE_OPTIONS:-DMaCh3_PYTHON_ENABLED=ON}

ARG INSTALL_OPTIONS
ENV INSTALL_OPTIONS=${INSTALL_OPTIONS:-"VERBOSE=1"}

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
RUN make ${INSTALL_OPTIONS} && make install

# KS: Need to set them here, otherwise container using this container will not be able to find MaCh3
ENV MaCh3_ROOT=${MACH3_INSTALL_DIR}
ENV PATH=${MaCh3_ROOT}/bin:${PATH} \
    LD_LIBRARY_PATH=${MaCh3_ROOT}/lib:${LD_LIBRARY_PATH}

# pip install pyMaCh3
WORKDIR ${MACH3_WORK_DIR}

#KS: If GPU is not being compiled this must mean ill configured cmake somewhere in MaCh3
RUN if ! mach3-config --has-feature GPU; then \
        echo "ERROR: MaCh3 was built without GPU support." >&2; \
        exit 1; \
    fi

RUN pip install .

# Start from MaCh3 install dir
WORKDIR ${MACH3_INSTALL_DIR}
