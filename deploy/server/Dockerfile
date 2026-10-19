# PDAC Tissue Explorer - app image for a self-hosted server (used by compose.yaml in this folder).
# Works on x86_64 and ARM64 servers: R and all packages come as prebuilt Ubuntu 24.04 packages
# (R 4.3.3, shiny 1.8, ggplot2 3.4, plotly 4.10), so nothing is compiled during the build.
# Build context is the repository root. Data files are not baked in; they live in a Docker volume
# filled once with get_data.R (see README.md in this folder).
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8 LC_ALL=C.UTF-8 TZ=Etc/UTC
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
       r-base-core r-cran-shiny r-cran-bslib r-cran-ggplot2 r-cran-plotly r-cran-dt \
       r-cran-matrix r-cran-markdown r-cran-patchwork r-cran-scales \
    && rm -rf /var/lib/apt/lists/*
# unprivileged user for the app (the base image already uses UID 1000 for "ubuntu")
RUN useradd --create-home shiny

WORKDIR /srv/pdac
COPY get_data.R /srv/pdac/get_data.R
COPY app/ /srv/pdac/app/
RUN mkdir -p /srv/pdac/app/data && chown -R shiny:shiny /srv/pdac

USER shiny
EXPOSE 3838
CMD ["R", "--quiet", "--no-save", "-e", "shiny::runApp('/srv/pdac/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
