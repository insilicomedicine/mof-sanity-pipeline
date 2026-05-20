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

# PLATON requires user-supplied license-restricted sources in code/. To build
# without PLATON entirely, pass `--build-arg WITHOUT_PLATON=1` — the runtime
# will detect the missing binary and skip PLATON in every code path.
ARG WITHOUT_PLATON=
COPY code/ /opt/project/code/
WORKDIR /opt/project/code
RUN if [ -n "$WITHOUT_PLATON" ]; then \
        echo "WITHOUT_PLATON=$WITHOUT_PLATON: skipping PLATON compilation"; \
        rm -f platon.f.gz xdrvr.c.gz platon.f xdrvr.c platon; \
    elif [ ! -f platon.f.gz ] || [ ! -f xdrvr.c.gz ]; then \
        echo "ERROR: PLATON sources not found in code/"; \
        echo "Comply with the PLATON license and download platon.f.gz and xdrvr.c.gz from https://www.platonsoft.nl/platon/pl030000.html"; \
        echo "and place them in the code/ directory before building."; \
        echo ""; \
        echo "To build without PLATON instead, pass --build-arg WITHOUT_PLATON=1"; \
        exit 1; \
    else \
        gzip -d platon.f.gz \
        && gzip -d xdrvr.c.gz \
        && gfortran -O3 -march=native -mtune=native --verbose \
                -static-libgfortran -static-libgcc \
                -Wl,-rpath,/opt/conda/lib \
                -o platon platon.f xdrvr.c \
                -L/usr/X11R6/lib -lX11 \
        && rm -f platon.f xdrvr.c; \
    fi


# =============================================================================
# Runtime image.
# (`final`) is declared at the end of the file as the default build target.
# =============================================================================
FROM sanity_runtime AS runtime

LABEL org.opencontainers.image.title="mof-sanity-pipeline" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="1.0.1"

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
# Optional smoke-test stage. Run with:
#     docker build --target smoke-test .
# Test fixtures are confined to this stage so the runtime image stays lean.
# =============================================================================
FROM runtime AS smoke-test
COPY test/test_cifs/ /opt/test_cifs/
# Per-CIF JSON fixtures live in /opt/test_fixtures and are shared between the
# two build modes (PLATON-only stage JSON gets skipped when WITHOUT_PLATON).
# CSV fixtures differ between modes — the no-platon variants are kept in a
# separate tree so the smoke script can pick the right one at runtime.
COPY test/fixtures/smoke/           /opt/test_fixtures/
COPY test/fixtures/smoke_no_platon/ /opt/test_fixtures_no_platon/

# Tunable smoke-test parameters. Override at build time, e.g.:
#   docker build --target smoke-test --build-arg SMOKE_NJOBS=8 .
#   docker build --target smoke-test --build-arg SMOKE_TOTAL_TIMEOUT=240 .
#   docker build --target smoke-test --build-arg WITHOUT_PLATON=1 .
ARG SMOKE_NJOBS=4
ARG SMOKE_TOTAL_TIMEOUT=600
ARG WITHOUT_PLATON=

# Each smoke step is its own RUN so BuildKit shows pass/fail per stage in the
# build progress output. State persists in /var/smoke between RUNs.
ENV DATA_DIR=/opt/test_cifs \
    FIXTURES_DIR=/opt/test_fixtures \
    FIXTURES_NO_PLATON_DIR=/opt/test_fixtures_no_platon \
    STATE_DIR=/var/smoke \
    SMOKE_NJOBS=${SMOKE_NJOBS} \
    SMOKE_TOTAL_TIMEOUT=${SMOKE_TOTAL_TIMEOUT} \
    WITHOUT_PLATON=${WITHOUT_PLATON}
RUN ln -sf /usr/local/bin/run_cli_smoke.sh /smoke

RUN /smoke prepare
RUN /smoke help
RUN /smoke check-help-content
RUN /smoke negative-unknown-flag
RUN /smoke negative-missing-input
RUN /smoke sanity-only
RUN /smoke check-stage-jsons
RUN /smoke oxichecker-only
RUN /smoke check-oxichecker-csv
RUN /smoke postprocess-only
RUN /smoke check-sanity-csv
RUN /smoke stage1-only
RUN /smoke single-file-input
RUN /smoke recursive
RUN /smoke custom-post-output
RUN /smoke no-obabel-run
RUN /smoke timeout-soft
RUN /smoke idempotency
RUN echo "ALL OK"

# =============================================================================
# Default target. Identical to `runtime`; placed last so `docker build .`
# without `--target` produces the runtime image, not the smoke-test one.
# =============================================================================
FROM runtime AS final
