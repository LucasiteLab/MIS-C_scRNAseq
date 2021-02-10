Kenneth B. Hoehn
2/1/2021
Scripts for performing BCR sequence analysis

Requirements:
Immcantation 4.0.0 Docker image
Change-O v1.0.1

R packages:
alakazam v1.1.0.999
shazam v1.0.2.999
dplyr v1.0.2
tidyr v1.1.2
stringr v1.4.0
stringdist v0.9.6.3
ggplot2 v3.3.2 
ggpubr v0.4.0

To reproduce analyses:
-Run the Immcantation docker image
-Run processData.sh to align BCR sequences to IMGT reference and convert to AIRR format.
-Run combineData.R to compile data into a single table with appropriate metadata
-Run cloneGermline.R to group BCRs into clonal clusters and reconstruct clonal germlines
-Run analysis.R to perform analyses and generate figures.
