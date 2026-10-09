# Optional: run the PDAC Tissue Explorer in a container (self-hosting, institutional servers, or any Docker host).
# Build from the repository root after placing the prepared data files in app/data/:
#   docker build -f deploy/docker/Dockerfile -t pdac-tissue-explorer .
#   docker run -p 3838:3838 pdac-tissue-explorer        # then open http://localhost:3838
# rocker/shiny pins CRAN to a dated Posit Package Manager snapshot, so package versions are reproducible.
FROM rocker/shiny:4.4.1

RUN install2.r --error --skipinstalled shiny bslib ggplot2 plotly DT Matrix markdown patchwork scales \
    && rm -rf /tmp/downloaded_packages

WORKDIR /srv/pdac-tissue-explorer
COPY app/ /srv/pdac-tissue-explorer/
USER shiny

EXPOSE 3838
CMD ["R", "--quiet", "-e", "shiny::runApp('/srv/pdac-tissue-explorer', host = '0.0.0.0', port = 3838)"]
