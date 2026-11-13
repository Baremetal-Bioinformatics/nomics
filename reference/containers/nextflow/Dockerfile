# Nextflow for linux/arm64. nextflow/nextflow publishes linux/amd64 only, so on
# arm64 `nomics verify fetch` builds this instead: the pinned image's own launcher
# and entry script (copied, not run), on the same JDK and OS it ships
# (Amazon Corretto 21.0.11, Amazon Linux 2023), with a static Docker CLI.
# NXF_IMAGE is the descriptor's pinned Nextflow image.
ARG NXF_IMAGE=docker.io/nextflow/nextflow:26.04.6
FROM --platform=linux/amd64 ${NXF_IMAGE} AS upstream
# Pinned by digest: this binary is handed the Docker socket.
FROM docker.io/library/docker:28.5.1-cli@sha256:9190b0613792e658a7783cf14b2d5ace5941bb68ede7276922ea36ee457d76ad AS dockercli
FROM docker.io/library/amazoncorretto:21.0.11-al2023
RUN yum install -y procps-ng shadow-utils which && yum clean all
ENV LANG=C.UTF-8 NXF_HOME=/.nextflow
COPY --from=upstream /usr/local/bin/nextflow /usr/local/bin/entry.sh /usr/local/bin/
COPY --from=dockercli /usr/local/bin/docker /usr/local/bin/docker
RUN mkdir /.nextflow && chmod 755 /usr/local/bin/nextflow /usr/local/bin/entry.sh && nextflow info
ENTRYPOINT ["/usr/local/bin/entry.sh"]
