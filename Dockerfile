FROM rocker/rstudio:latest
MAINTAINER Joerg Lenhard <joerg.lenhard@kau.se>

COPY /data /home/rstudio/arch-metrics-replication/data
COPY /scripts /home/rstudio/arch-metrics-replication/scripts