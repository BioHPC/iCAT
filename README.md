## iCAT: Immune Cells Analysis software Tool

High throughput sequencing of immune cell receptor sequences presents a unique opportunity to inform our understanding of immunological response to infection and how we can detect it. However, the nature of this process requires analysis of massive datasets that requires streamlining. Here we present ICAT, an R-based method for tackling this data and applying statistical analyses to find significant patterns and divergences.

<br/>

![Alt text](/screenshot/icat.png?raw=true "Screeshot")

<br/>

### Pre-requisites:

* R version
  * Download R (>3.4.0) version from CRAN.
    * Windows: https://cran.r-project.org/bin/windows/base/
    * Mac OS X: https://cran.r-project.org/bin/macosx/
    * Linux: https://cran.r-project.org/bin/linux/

* Libraries:
    - devtools

To install devtools, use the command:

        > install.packages("devtools") 

### Download/Install iCAT Package:

        > devtools::install_github("BioHPC/iCAT") 


### Run iCAT Package:

        > library(iCAT)
        > iCATinteractive()

