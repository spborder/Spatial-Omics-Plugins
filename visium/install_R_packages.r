#!/usr/bin/env Rscript
library(stringr)
## This script is used for installing a series of R packages and returning an error if a package fails to install.
## (will also trigger an error if a dependency is not installed)
packages <- commandArgs(trailingOnly=TRUE)
package_type <- packages[1]
packages <- packages[2:length(packages)] 

if (package_type=="install-packages") {
    for (l in packages) {
        install.packages(l, dependencies=TRUE, repos = "https://cran.rstudio.com/");

        if ( ! library(l, character.only=TRUE,logical.return=TRUE) ) {
            quit(status=1, save='no')
        }

    }
} else if (package_type=='BiocManager') {
    for (l in packages) {
        BiocManager::install(l, dependencies=TRUE,update=TRUE,ask=FALSE);

        if ( ! library(l, character.only=TRUE,logical.return=TRUE) ) {
            quit(status=1, save='no')
        }

    }
} else if (package_type=='GitHub') {

    for (l in packages) {
        split_package_name <- str_split(l,"&")
        remotes::install_github(l[1],ref=l[2], dependencies=TRUE)

        if ( ! library(l, character.only=TRUE,logical.return=TRUE) ) {
            quit(status=1, save='no')
        }

    }

}

