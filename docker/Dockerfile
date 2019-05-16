# Galaxy - PIPE_T
#
# VERSION       0.2

FROM bgruening/galaxy-stable:latest

MAINTAINER D. CANGELOSI, davidecangelosi@gaslini.org

ENV GALAXY_CONFIG_BRAND="Galaxy PIPE-T"
ENV ENABLE_TTS_INSTALL True
ENV CONDA_ENSURE_CHANNELS="r,iuc,conda-forge,bioconda,defaults"
ENV DEPENDENCY_RESOLVERS_CONFIG_FILE="config/dependency_resolvers_conf.xml"


RUN apt-get -y update
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev curl libcairo2-dev

#Workaround for ubuntu
#RUN cp /etc/ssl/private/ssl-cert-snakeoil.key /etc/
#RUN chown root:ssl-cert /etc/ssl-cert-snakeoil.key
#RUN sed -i -e "s/\/ssl\/private//g" /etc/postgresql/9.3/main/postgresql.conf
ADD docker/dependency_resolvers_conf.xml $GALAXY_ROOT/config/dependency_resolvers_conf.xml

#install pipe_t from toolshed
ADD docker/tools.yml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs --yes && \  
    rm /export/galaxy-central/ -rf 

#RUN chown galaxy:galaxy /usr/local/lib/R/site-library

# Mark folders as imported from the host.
#VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]
