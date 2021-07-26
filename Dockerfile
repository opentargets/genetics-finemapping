FROM ubuntu:21.04

ENV DEBIAN_FRONTEND=noninteractive

# Create non-root user 
ARG UID
ARG GID
RUN groupadd -g $GID -o otg
RUN useradd -m -u $UID -g $GID -o -s /bin/bash otg

# Do everything that requires root user
# Install dependencies
RUN apt-get update && \
    apt-get install --yes ant openjdk-8-jdk wget bzip2 parallel unzip && \
    apt-get clean

# fix certificate issues (as root)
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# install micromamba
RUN mkdir -p /software/micromamba && \
    cd /software/micromamba && \
    wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/0.15.2 | tar -xvj bin/micromamba && \
    chown -R otg:otg /software/micromamba
ENV PATH="/software/micromamba/bin:${PATH}"

# download gcta
RUN wget https://cnsgenomics.com/software/gcta/bin/gcta_1.92.3beta3.zip -P /software/gcta && \
    unzip /software/gcta/gcta_1.92.3beta3.zip -d /software/gcta && \
    rm /software/gcta/gcta_1.92.3beta3.zip && \
    chown -R otg:otg /software/gcta
ENV PATH="/software/gcta/gcta_1.92.3beta3:${PATH}"

# switch to otg user
USER otg

# create conda/mamba environment
COPY ./environment.yaml /home/otg
RUN micromamba install --name base --file /home/otg/environment.yaml --root-prefix /software/micromamba --yes

# set JAVA_HOME (useful for docker commandline)
ENV JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/"

# install Google Cloud SDK
# RUN curl https://sdk.cloud.google.com | bash

# copy all files of the repo and change owner
COPY ./ /finemapping
USER root
RUN chown -R otg:otg /finemapping
USER otg

# set default directory
WORKDIR /finemapping

# default command
CMD ["/bin/bash"]
