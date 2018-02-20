FROM rocker/rstudio:3.4.3
MAINTAINER Joerg Lenhard <joerg.lenhard@kau.se>

COPY /data /home/rstudio/arch-metrics-replication/data
COPY /scripts /home/rstudio/arch-metrics-replication/scripts