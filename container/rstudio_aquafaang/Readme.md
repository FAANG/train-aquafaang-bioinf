# R Studio with all libraries for the AquaFaang course pre-installed

Links to sources:
* [Bioconductor](http://bioconductor.org/)
* [ChIPseeker](https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
* [DiffBind](https://www.bioconductor.org/packages/release/bioc/html/DiffBind.html)
* [GenomicFeature](https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html])
* [Tidyverse](https://www.tidyverse.org/)

## Run
```bash
 docker run \
     -e PASSWORD=bioc \
     -p 8787:8787 \
     --mount type=bind,source="$(pwd)",target=/mnt \
     juettemann/rstudio_aquafaang
```

Open Browser
* [http://localhost:8787](http://localhost:8787)
* user: rstudio
* password: bioc (unless changed in the docker command)

By default, files in the directory where the docker command is run are
accessible in /mnt

## Installing dependencies without docker

* [R/studio](https://www.rstudio.com/products/rstudio/download/)
* [Bioconductor](https://www.bioconductor.org/install/)

Once Rstudio and Bioconductor are installed:

```r
install.packages("tidyverse")
BiocManager::install("DiffBind")
```

For Ubuntu users, some system dependencies must also be installed:

```bash
sudo apt install librsvg2-dev libv8-dev libcurl4-openssl-dev libssl-dev libxml2-dev
```

