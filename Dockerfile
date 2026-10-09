FROM mambaorg/micromamba:2.9.0

# Keep the same environment as the native installation.
ENV ENV_NAME=cembio_eurac \
    RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_CACHE_SYMLINKS=FALSE
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create --yes --file /tmp/environment.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
RUN mkdir -p /workspace /opt/renv /opt/cembio /cache && \
    chown -R "$MAMBA_USER:$MAMBA_USER" /workspace /opt/renv /opt/cembio && \
    chmod 1777 /cache
USER $MAMBA_USER
WORKDIR /workspace

# Restore before copying workflow code so edits reuse the dependency layer.
COPY --chown=$MAMBA_USER:$MAMBA_USER DESCRIPTION environment.yml renv.lock .Rprofile ./
COPY --chown=$MAMBA_USER:$MAMBA_USER renv/activate.R renv/activate.R
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/bootstrap_environment.R scripts/bootstrap_environment.R
ARG MAKEFLAGS=-j2
RUN Rscript scripts/bootstrap_environment.R && \
    Rscript -e 'writeLines(paste0(".libPaths(", encodeString(renv::paths$library(), quote = "\""), ")"), "/opt/cembio/Rprofile")' && \
    chmod -R a+rX /opt/renv && \
    rm -rf /opt/renv/cache

# Runtime libraries live outside the bind-mounted checkout. The fixed profile
# uses those libraries without activating the host's project-local renv.
ENV R_PROFILE_USER=/opt/cembio/Rprofile \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    HOME=/tmp \
    XDG_CACHE_HOME=/cache \
    BFC_CACHE=/cache/BiocFileCache
COPY --chown=$MAMBA_USER:$MAMBA_USER . .
COPY docker/run-workflow.sh /usr/local/bin/run-workflow
USER root
RUN chmod 755 /usr/local/bin/run-workflow
USER $MAMBA_USER
RUN Rscript scripts/check_environment.R && \
    XDG_CACHE_HOME=/tmp/cembio-smoke-cache \
    quarto render docker/smoke-test.qmd --output-dir /tmp/cembio-smoke && \
    rm -rf /tmp/cembio-smoke /tmp/cembio-smoke-cache docker/smoke-test_files

# Permit writes to mounted checkouts by default; Linux users can supply --user
# to create results with their own UID/GID.
USER root
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/usr/local/bin/run-workflow"]
CMD ["help"]
