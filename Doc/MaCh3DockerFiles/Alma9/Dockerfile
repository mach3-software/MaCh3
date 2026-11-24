# To run use:  docker build -t mach3 .
# Stage 1: build MaCh3, bu having separate stage we can reduce hosting space
# KS: Get glorious container from Kamil which will work as a base
FROM kamilskwarczynski/nukamil:alma9root64004 AS mach3_build

ARG MACH3_VERSION
ENV MACH3_VERSION=${MACH3_VERSION:-develop}

ARG CMAKE_OPTIONS
ENV CMAKE_OPTIONS=${CMAKE_OPTIONS:-DMaCh3_PYTHON_ENABLED=ON}

ARG INSTALL_OPTIONS
ENV INSTALL_OPTIONS=${INSTALL_OPTIONS:-"VERBOSE=1"}

ARG REPOSITORY_URL
ENV REPOSITORY_URL=${REPOSITORY_URL:-https://github.com/mach3-software/MaCh3}

ENV MACH3_SRC_DIR=/tmp/MaCh3-src
ENV MACH3_BUILD_DIR=/tmp/MaCh3-build
ENV MACH3_INSTALL_DIR=/opt/MaCh3

# Clone, build out-of-tree, install to a clean prefix
RUN git clone ${REPOSITORY_URL} ${MACH3_SRC_DIR} \
 && git -C ${MACH3_SRC_DIR} checkout ${MACH3_VERSION} \
 && cmake -S ${MACH3_SRC_DIR} -B ${MACH3_BUILD_DIR} \
          -DCMAKE_INSTALL_PREFIX=${MACH3_INSTALL_DIR} \
          ${CMAKE_OPTIONS} \
 && make -C ${MACH3_BUILD_DIR} -j"$(nproc)" ${INSTALL_OPTIONS} \
 && make -C ${MACH3_BUILD_DIR} install \
 && find ${MACH3_INSTALL_DIR} -type f -name '*.so*' -exec strip --strip-unneeded {} + || true

# Build the python bindings as a wheel (needs source + installed MaCh3)
ENV MaCh3_ROOT=${MACH3_INSTALL_DIR} \
    PATH=${MACH3_INSTALL_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH=${MACH3_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}
RUN pip wheel --no-cache-dir -w /wheels ${MACH3_SRC_DIR}

# Stage 2: this is actual image we gonna upload,
# Add a label for the author
LABEL maintainer="The MaCh3 Collaboration"
LABEL website="https://mach3-software.github.io/MaCh3/"
LABEL compiler="GNU 11.5.0"
LABEL root_version="v6.40.04"
LABEL python_version="v3.12.16"
LABEL org.opencontainers.image.description="Official MaCh3 container"

# Declare argument
ARG MACH3_VERSION
ENV MACH3_VERSION=${MACH3_VERSION:-develop}

ENV MACH3_INSTALL_DIR=/opt/MaCh3
COPY --from=mach3_build /opt/MaCh3 /opt/MaCh3

# KS: Need to set them here, otherwise containers built on this one can't find MaCh3
ENV MaCh3_ROOT=${MACH3_INSTALL_DIR} \
    PATH=${MACH3_INSTALL_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH=${MACH3_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}

# Install the wheel
RUN --mount=type=bind,from=mach3_build,source=/wheels,target=/wheels \
    pip install --no-cache-dir /wheels/*.whl

WORKDIR ${MACH3_INSTALL_DIR}
