#docker build -t alma9root64004 .
#docker tag alma9root64004:latest kamilskwarczynski/nukamil:alma9root64004
#docker push kamilskwarczynski/nukamil:alma9root64004

########################################
# Stage 0: Build dependencies (devel packages, NOT shipped)
########################################
FROM almalinux:9 AS build-deps

RUN dnf -y install --setopt=install_weak_deps=False --nodocs epel-release && \
    dnf -y --enablerepo=crb install --setopt=install_weak_deps=False --nodocs \
        wget git gl2ps gl2ps-devel \
        gcc gcc-c++ gcc-gfortran \
        make cmake openssl-devel \
        perl giflib-devel xrootd-devel xrootd-client-devel \
        which findutils procps-ng libSM-devel libICE-devel \
        ncurses ncurses-compat-libs libuuid-devel \
        tar gzip xz unzip krb5-devel fontconfig-devel \
        python3.12 python3.12-devel python3.12-pip python3.12-numpy \
        libarrow-devel \
        openssh-clients fftw-devel gsl-devel gmp-devel hdf5-devel cfitsio-devel graphviz-devel \
        libX11-devel libXpm-devel libXft-devel libXext-devel libXt-devel libXaw-devel \
        mesa-libGLU-devel glew-devel motif-devel \
        libpng-devel libjpeg-turbo-devel freetype-devel \
        sqlite-devel libxml2-devel davix-devel \
        lz4-devel zlib-devel xz-devel pcre-devel \
        tbb-devel protobuf-devel protobuf-compiler && \
    dnf clean all && rm -rf /var/cache/dnf

########################################
# Stage 1: Build ROOT against Python 3.12
########################################
FROM build-deps AS root-builder

ARG ROOT_VERSION=v6-40-04
ARG NCORES=16

RUN git clone --depth 1 --branch ${ROOT_VERSION} \
        https://github.com/root-project/root.git /tmp/root-src && \
    mkdir /tmp/root-build

WORKDIR /tmp/root-build

RUN cmake /tmp/root-src \
        -DCMAKE_INSTALL_PREFIX=/opt/root/${ROOT_VERSION} \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -Dminuit2_omp=OFF \
        -Darrow=ON \
        -Dmathmore=ON \
        -Dclad=OFF \
        -Dbuiltin_zstd=ON \
        -Droot7=OFF \
        -DPython3_EXECUTABLE=/usr/bin/python3.12 && \
    make -j${NCORES} && make install

# Slim the install: drop tutorials, strip symbols
RUN rm -rf /opt/root/${ROOT_VERSION}/tutorials && \
    find /opt/root/${ROOT_VERSION} -type f \
        \( -name '*.so' -o -name '*.so.*' -o -path '*/bin/*' \) \
        -exec strip --strip-unneeded {} + 2>/dev/null || true

########################################
# Stage 2: Final runtime image (runtime libs only, no -devel)
########################################
FROM almalinux:9 AS final

ARG ROOT_VERSION=v6-40-04
# Python packages from pip (drop pyarrow if you don't need it, it is large)
ARG PIP_PACKAGES="pyyaml matplotlib pyarrow"

RUN dnf -y install --setopt=install_weak_deps=False --nodocs epel-release && \
    dnf -y --enablerepo=crb install --setopt=install_weak_deps=False --nodocs \
        git \
        gcc gcc-c++ gcc-gfortran make cmake \
        which findutils procps-ng less nano vim-minimal tmux openssh-clients tar gzip xz unzip \
        python3.12 python3.12-devel python3.12-pip python3.12-numpy \
        gl2ps giflib xrootd-libs xrootd-client-libs \
        libSM libICE libuuid ncurses-compat-libs krb5-libs fontconfig \
        libX11 libXpm libXft libXext libXt libXaw \
        mesa-libGLU libglvnd-opengl glew motif libpng libjpeg-turbo freetype \
        sqlite-libs libxml2 davix-libs lz4-libs zlib xz-libs pcre \
        tbb tbb-devel protobuf libarrow \
        fftw-libs gsl gmp hdf5 cfitsio graphviz && \
    python3.12 -m pip install --no-cache-dir ${PIP_PACKAGES} && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python3 && \
    ln -sf /usr/bin/python3.12 /usr/local/bin/python && \
    ln -sf /usr/bin/pip3.12   /usr/local/bin/pip3 && \
    ln -sf /usr/bin/pip3.12   /usr/local/bin/pip && \
    dnf clean all && rm -rf /var/cache/dnf /root/.cache

# copy ROOT
COPY --from=root-builder /opt/root/${ROOT_VERSION} /opt/root/${ROOT_VERSION}

# environment
ENV SHELL=/bin/bash
ENV ROOTSYS=/opt/root/${ROOT_VERSION}
ENV PATH=/opt/root/${ROOT_VERSION}/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/root/${ROOT_VERSION}/lib:/opt/root/${ROOT_VERSION}/lib/root
ENV PYTHONPATH=/opt/root/${ROOT_VERSION}/lib/root:/opt/root/${ROOT_VERSION}/lib

SHELL ["/bin/bash", "-c"]

# interactive shells get the ROOT environment automatically
RUN echo "source ${ROOTSYS}/bin/thisroot.sh" >> /etc/bashrc

# verify toolchain and catch any missing runtime libs
RUN source ${ROOTSYS}/bin/thisroot.sh && \
    root-config --version && \
    python3 --version && \
    pip --version && \
    python3 -c "import ROOT; print(ROOT.gROOT.GetVersion())" && \
    ! ( find ${ROOTSYS}/lib -name '*.so*' -exec ldd {} + 2>/dev/null | grep "not found" | sort -u | grep . )

CMD ["/bin/bash"]
