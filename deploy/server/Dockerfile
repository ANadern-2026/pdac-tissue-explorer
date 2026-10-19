# PDAC Tissue Explorer - app image for a self-hosted server (used by compose.yaml in this folder).
# Build context is the repository root. Data files are not baked in; they live in a Docker volume
# filled once with get_data.R (see README.md in this folder).
FROM rocker/shiny:4.4.1

RUN install2.r --error --skipinstalled shiny bslib ggplot2 plotly DT Matrix markdown patchwork scales \
    && rm -rf /tmp/downloaded_packages

WORKDIR /srv/pdac
COPY get_data.R /srv/pdac/get_data.R
COPY app/ /srv/pdac/app/
RUN mkdir -p /srv/pdac/app/data && chown -R shiny:shiny /srv/pdac

USER shiny
EXPOSE 3838
CMD ["R", "--quiet", "-e", "shiny::runApp('/srv/pdac/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
