FROM bioconductor/bioconductor_docker:RELEASE_3_12

RUN Rscript -e "install.packages(c('rlang', 'knitr' ), repos = 'https://cloud.r-project.org')"

RUN Rscript -e "BiocManager::install()"
RUN Rscript -e "BiocManager::install('DiffBind')"
RUN Rscript -e "BiocManager::install('ChIPseeker')"
RUN Rscript -e "BiocManager::install('GenomicFeatures')"
RUN Rscript -e "BiocManager::install('tidyverse')"
