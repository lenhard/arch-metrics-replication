# Replication package for the paper: "Exploring the Suitability of Source Code Metrics for Indicating Architectural Inconsistencies"

The following folders and files are part of this replication package:
 - `data\`: Contains the raw metric data used in the study in the form of comma-separated-value files. There is one file for each case study system, e.g. `jabref-data.csv`.
 - `models\`: Contains textual version and images of the modules of the systems and the mapping of the source code to the indended module structure, e.g. `jabref-archmodel.txt` and `jabref-model.png`.
 - `scripts\`: Contains R scripts that have been used to execute the statistical computations discussed in the paper. There is one file for each case study system, e.g. `jabref-analysis.R`.
 
To reproduce the calculations, simply start your favourite R development environment and open the R scripts. The scripts assume that the repository is located on your R working path. More instructions can be found inside the scripts.

For convenience, the repository contains a Dockerfile which sets up RStudio in a docker container. The advantage of that is that you do not have to install the environment locally, but can access everthing through your brower. To run this repository inside a container, you need to have Docker installed on your system. If you do, `cd` to the directory of the repository and do the following:
 - build an image: `docker build -t arch-metrics-replication .`
 - run the image (RStudio will be available on port 8787): `docker run --rm -p 8787:8787 arch-metrics-replication`
 - open your browser at http://localhost:8787/ (Linux) or http://YOUR-DOCKER-MACHINE-IP:8787/
 - RStudio's default credentials are username: `rstudio`, password: `rstudio`
 - Open one of the R scripts in folder `arch-metrics-replication`, read the code and execute it


