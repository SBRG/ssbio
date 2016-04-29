FROM andrewosh/binder-base

USER root

# Add R dependencies
RUN apt-get update
RUN apt-get install -y r-base libzmq3-dev

COPY install-irkernel.R /home/install-irkernel.R

RUN R --no-save < /home/install-irkernel.R

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")'

USER main