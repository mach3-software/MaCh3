#To run use:  docker build -t mach3 .
#FROM rootproject/root:6.24.06-centos7 AS mach3_build
#KS: Get glorious container from Luke which will work as a base
FROM picker24/root_v6_26_10:alma9 AS mach3_build

ENV MACH3_VERSION="develop"
ENV MACH3_INSTALL_DIR=/opt/mach3/${MACH3_VERSION}

RUN mkdir -p /opt/MaCh3/

WORKDIR /opt/
RUN --mount=type=ssh git clone https://github.com/mach3-software/MaCh3 mach3-src
WORKDIR /opt/mach3-src
RUN git checkout ${MACH3_VERSION}

RUN mkdir -p ${MACH3_INSTALL_DIR}
WORKDIR ${MACH3_INSTALL_DIR}
RUN cmake /opt/mach3-src
RUN make VERBOSE=1 && make install

RUN source ${MACH3_INSTALL_DIR}/bin/setup.MaCh3.sh

WORKDIR ${MACH3_INSTALL_DIR}
