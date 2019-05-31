# Immune Cells Analysis software Tool <img src="man/figures/logo.png" align="right" width="135"/> 

High throughput sequencing of immune cell receptor sequences presents a unique opportunity to inform our understanding of immunological response to infection and how we can detect it. However, the nature of this process requires analysis of massive datasets that requires streamlining. Here we present ICAT, an R-based method for tackling this data and applying statistical analyses to find significant patterns and divergences.

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
        
*Note*: iCAT also uses shiny, shinyjs, data.table, ggplot2, DT, hash, and magrittr. However, those packages will be installed if using `install_github` from below.

### Installing iCAT

Using an R interface, type:
```
> devtools::install_github("BioHPC/iCAT") 
```
