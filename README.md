# iCAT: Immune Cells Analysis software Tool

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
        
*Note*: iCAT also uses shiny, shinyjs, data.table, ggplot2, DT, hash, and magrittr. However, those packages will be installed if using `install_github` from below.

## Installing iCAT

Using an R interface, type:

        > devtools::install_github("BioHPC/iCAT") 


## Graphical Interface Workflow

In R:

        > library(iCAT)
        > iCATinteractive()

This will launch a graphical user interface (GUI) for iCAT. The GUI has three tabs, separating major functionalities: training, library, and prediction.

**Training:**

1) In the _Training_ tab, enter your negative training samples (naïve, unexposed, uninfected, pre-infection, etc) using the `Browse` button.

Individual samples’ sequencing data should be in .tsv format.

2) Repeat step 1 for positive training samples (exposed, infected, etc.)

3) Choose if you want to analyze data by: 
- `CDR3 Amino Acid Sequence` (TCRs will need the same CDR3 region to be called ‘Identical’)
- `TCRV-CDR3-TCRJ` (TCRs will need the same TCRBV segment, CDR3 region, and TCRJ segment to be called ‘Identical’) *Recommended*
- `Nucleic Acid (DNA)` (TCRs will need the exact same DNA rearrangements/sequence across TCRBV, CDR3, and TCRJ)

4) Choose the `Max p-value`, which determines the minimal degree of statistical significance that iCAT will accept as being potentially "associated" with the positive group. Defaults to _p_ < 0.1.

5)	Choose the `Min Threshold of Public Sequences`, which determines the minimum number of training samples a TCR sequence must be observed in to be considered as potentially "associated" with the positive group. Defaults to 1. Recommend setting at 10% of positive training samples. E.g. if there are 30 positive training samples, a recommended minimum threshold is 3. 

6) Once all options are selected click `Train Model`



