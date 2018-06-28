# iCAT: Immune Cells Analysis software Tool

PI: Dr. Richard DiPaolo
PI: Dr. Tae-Hyuk (Ted) Ahn

High throughput sequencing of immune cell receptor sequences presents a unique opportunity to inform our understanding of immunological response to infection and how we can detect it. However, the nature of this process requires analysis of massive datasets that requires streamlining. Here we present ICAT, an R-based method for tackling this data and applying statistical analyses to find significant patterns and divergences.

# Pre-requisites:

* R version
  * Download R (>3.4.0) version from CRAN.
    * Windows: https://cran.r-project.org/bin/windows/base/
    * Mac OS X: https://cran.r-project.org/bin/macosx/
    * Linux: https://cran.r-project.org/bin/linux/

- data.table package

  * install by using the following R command:

        > install.packages("data.table")  

- dplyr package

  * install by using the following R command:

        > install.packages("dplyr")  

- plyr package

  * install by using the following R command:

        > install.packages("plyr")  

# Download iCAT Package:

* Click the "Clone or download" button located at top-right
* Click "Download ZIP" and extract the compressed ZIP file
* Or you can clone with Git (https://help.github.com/articles/cloning-a-repository/)

# Installing iCAT Package:

* You can run the R script without installation 

# Run iCAT Package:

* Provide pre-vaccination tab-delimited files into Pre directory
* Provide post-vaccination tab-delimited files into Post directory
* Run R script as below
* source("iCAT.r")

