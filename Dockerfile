# syntax=docker/dockerfile:1.7

# =============================================================================
# Runtime base: conda + pinned Python stack + openbabel runtime + gosu.
# Used as the base of the final runtime image.
# =============================================================================
FROM continuumio/miniconda3:25.1.1-2 AS sanity_runtime

# Allow BuildKit cache mounts to retain apt downloads between builds.
RUN rm -f /etc/apt/apt.conf.d/docker-clean

# Runtime apt deps only.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
        openbabel \
        gosu \
        ca-certificates

# Conda packages.
RUN conda install -y -n base -c conda-forge \
        pymatgen=2025.2.18 \
        backports=1.0 \
        backports.cached-property=1.0.2 \
    && conda clean -afy

# Pip packages pinned. libconeangle 0.1.2 has no prebuilt wheel for cp312/linux
# and is built from source via scikit-build, so we install build tools
# transiently inside the same RUN and purge them at the end to keep the
# runtime layer lean.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    --mount=type=cache,target=/root/.cache/pip \
    apt-get update \
    && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        gfortran \
        cmake \
        ninja-build \
    && pip install \
        structuregraph_helpers==0.0.9 \
        ase==3.24.0 \
        cached_property==2.0.1 \
        element_coder==0.0.8 \
        libconeangle==0.1.2 \
        openbabel-wheel==3.1.1.21 \
        rdkit==2024.9.5 \
    && apt-get purge -y --auto-remove \
        gcc \
        g++ \
        gfortran \
        cmake \
        ninja-build

# =============================================================================
# Build base: runtime + compilers, -dev headers, CMake.
# =============================================================================
FROM sanity_runtime AS sanity_build

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        gfortran \
        make \
        unzip \
        gzip \
        wget \
        libopenbabel-dev \
        libboost-filesystem-dev \
        libboost-system-dev \
        libboost-thread-dev \
        libeigen3-dev \
        libx11-dev

# CMake newer than the bookworm package; verified against the published SHA256.
ARG CMAKE_VERSION=3.27.0
ARG CMAKE_SHA256=b4cdcf94d06cd04e065cb0607535d76a4c12c167a6f99d3f4dae31f09bedb77c
RUN wget -q https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh \
    && echo "${CMAKE_SHA256}  cmake-${CMAKE_VERSION}-linux-x86_64.sh" | sha256sum -c - \
    && sh cmake-${CMAKE_VERSION}-linux-x86_64.sh --skip-license --prefix=/usr/local \
    && rm cmake-${CMAKE_VERSION}-linux-x86_64.sh

# =============================================================================
# Builder: compiles libcif and PLATON. Selective COPYs maximise cache hits.
# =============================================================================
FROM sanity_build AS builder

# libcif rebuilds only when libcif/ changes.
COPY libcif/ /opt/project/libcif/
WORKDIR /opt/project/libcif
RUN cmake -S . -B build/ --fresh -DLIBCIF_BLAKE2=ON \
    && cmake --build build/ -j8 \
    && mkdir -p /opt/libcif/bin \
    && cp build/examples/cif_to_building_blocks_v2 /opt/libcif/bin/ \
    && cp build/examples/cif_to_building_blocks_v3 /opt/libcif/bin/ \
    && cp build/examples/cif_check                 /opt/libcif/bin/

# PLATON requires user-supplied license-restricted sources in code/.
# Static gfortran/gcc linkage + rpath to conda's libquadmath so the binary
# runs in the runtime stage (which has gfortran purged for image size).
COPY code/ /opt/project/code/
WORKDIR /opt/project/code
RUN if [ ! -f platon.f.gz ] || [ ! -f xdrvr.c.gz ]; then \
        echo "ERROR: PLATON sources not found in code/"; \
        echo "Comply with the PLATON license and download platon.f.gz and xdrvr.c.gz from https://www.platonsoft.nl/platon/pl030000.html"; \
        echo "and place them in the code/ directory before building."; \
        exit 1; \
    fi \
    && gzip -d platon.f.gz \
    && gzip -d xdrvr.c.gz \
    && gfortran -O3 -march=native -mtune=native --verbose \
            -static-libgfortran -static-libgcc \
            -Wl,-rpath,/opt/conda/lib \
            -o platon platon.f xdrvr.c \
            -L/usr/X11R6/lib -lX11 \
    && rm -f platon.f xdrvr.c


# =============================================================================
# Runtime image.
# (`final`) is declared at the end of the file as the default build target.
# =============================================================================
FROM sanity_runtime AS runtime

LABEL org.opencontainers.image.title="mof-sanity-pipeline" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="paper"

# Build artefacts from the builder stage.
COPY --from=builder /opt/libcif/bin    /opt/libcif/bin
COPY --from=builder /opt/project/code/ /code/

# Project files copied straight from the build context — keeps libcif/PLATON
# build cache valid when only Python sources or scripts change.
COPY oxichecker/              /opt/oxichecker/
COPY generate_entrypoint.sh   /usr/local/bin/generate_entrypoint.sh
COPY scripts/run_cli_smoke.sh /usr/local/bin/run_cli_smoke.sh
COPY docker-entrypoint.sh     /usr/local/bin/docker-entrypoint.sh

RUN chmod +x \
        /usr/local/bin/generate_entrypoint.sh \
        /usr/local/bin/run_cli_smoke.sh \
        /usr/local/bin/docker-entrypoint.sh \
    && /usr/local/bin/generate_entrypoint.sh \
    && chmod 755 /entrypoint.sh

# Optional CLI smoke tests baked into the build. Pass:
#   docker build --build-arg RUN_ENTRYPOINT_TESTS=yes .
# Test CIFs are baked in only when this flag is enabled.
ARG RUN_ENTRYPOINT_TESTS=no
ENV RUN_ENTRYPOINT_TESTS=${RUN_ENTRYPOINT_TESTS}
COPY test/test_cifs/ /opt/test_cifs/
RUN set -e; \
    case "$RUN_ENTRYPOINT_TESTS" in \
      [Yy][Ee][Ss]|[Tt][Rr][Uu][Ee]|1) DATA_DIR=/opt/test_cifs /usr/local/bin/run_cli_smoke.sh;; \
      *) echo "Skipping CLI smoke tests (RUN_ENTRYPOINT_TESTS=$RUN_ENTRYPOINT_TESTS)";; \
    esac

# Pre-compile Python sources to bytecode for faster startup.
RUN python3 -m compileall /opt/oxichecker/ /code/ -b \
    && find /opt/oxichecker/ /code/ -name "*.py" -not -name "__*" \
            -exec python3 -O -m py_compile {} \;

ENV NUMEXPR_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    TLIMIT=10000000 \
    PYTHONPATH=/opt/oxichecker:/code \
    CIF_TO_BB_PATH=/opt/libcif/bin/cif_to_building_blocks_v3 \
    LIBCIF_BINARY_PATH=/opt/libcif/bin/cif_check \
    MPLCONFIGDIR=/var/cache/matplotlib

# Writable data + matplotlib cache directories for arbitrary host UIDs. Sticky
# bit so users can only remove their own files (same model as /tmp); replaces
# blanket chmod 777. MPLCONFIGDIR is set because the gosu'd runtimeuser has no
# $HOME, which otherwise triggers a matplotlib fallback warning at every import.
RUN mkdir -p /data /var/cache/matplotlib \
    && chmod 1777 /data /var/cache/matplotlib

# NOTE: container runs as root on purpose. docker-entrypoint.sh invokes `gosu`
# to drop privileges to a host-mapped user at startup; switching USER here
# would break that mechanism (gosu must be called from root).
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]


# =============================================================================
# Default target. Identical to `runtime`; placed last so `docker build .`
# without `--target` produces the runtime image.
# =============================================================================
FROM runtime AS final
