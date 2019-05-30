# iCAT: Immune Cells Analysis software Tool <img src="inst/app/www/cat2.png" align="right" width="125"/>


High throughput sequencing of immune cell receptor sequences presents a unique opportunity to inform our understanding of immunological response to infection and how we can detect it. However, the nature of this process requires analysis of massive datasets that requires streamlining. Here we present ICAT, an R-based method for tackling this data and applying statistical analyses to find significant patterns and divergences.

<br/>


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

![Alt text](/screenshot/icat.png?raw=true "Training")


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

A progress bar will show on the bottom-right corner to update on the satatus of training. After finishing, the _training_ tab will show some exploratory tables and a figure regarding the training data and the model built, which can all be downloaded to the user's machine. In addition, the _library_ and _prediction_ tab will unlock.

**Library:**

![Alt text](/screenshot/lib-icat.png?raw=true "Library")

The _Library_ tab displays a table consisting of the "target associated receptor sequences" (TARS), determined to be statistically associated with exposure to the target/agent/pathogen. The table displays each sequence, number of positive and negative training samples the sequence is present in/absent from, and how statistically associated the sequence is to the positive training data (*p*-value). The table can be downloaded to the user's computer for further analysis using excel, commandline, etc using the download `Table` button below the table.

**Prediction:**

![Alt text](/screenshot/pred-icat.png?raw=true "Prediction")

The _Prediction_ tab allows the user to add sequencing data from unknown samples (e.g. not included in the previous training data) for classification as "Positive" or "Negative" and determining the accuracy of the diagnostic assay.

1)	Use the `Browse` button to add samples for prediction. Multiple samples may be uploaded simultaneously.
2)	Click `Predict Independent Sample`.

A table will appear after analysis is complete. The table displays sample names along with the prediction "Positive" (red)/ "Negative" (blue), and displays the ‘%TARS’: the percent of individual sequences from the sample that are included in the TARS library. The prediction results can be downloaded as a table.

## R-interface Workflow

After loading iCAT with `library(iCAT)`

1) Define the parameters you would like to use for building the model:
```
> FIELD <- "vGeneName aminoAcid jGeneName"
> P_CUTOFF <- 0.1
> MIN_PUBLIC <- 2
```     
2) Make lists of .tsv Positive and Negative training samples:

```     
> listPos <- tsvDir("path/to/positve/samples/")
> listNeg <- tsvDir("path/to/negative/samples/")
       
```     
 - _optional_ Collect summary statistics about training samples:

```     
> trnStats(listPos, listNeg, FIELD)
#>         # Samples # Clonotypes # Unique Sequences
#>Negative         9       101942             281336
#>Positive         9        34158              89150
```     
3) Read in Positive and Negative training samples:

```     
> naive <- readTrn(listNeg, FIELD, "naive")
> vaccs <- readTrn(listPos, FIELD, "vacc")       
```     
4) Build a model using the training data:
```      
> mod <- train(naive, vaccs, listNeg, listPos, FIELD, P_CUTOFF, MIN_PUBLIC, NULL)
```     
       
 - _optional_ Produce a table estimating the classification accuracy of the model: 

```     
> classMat(mod)
```     
 - _optional_ Produce a figure showing % of TCR associated with positve samples in positive and negative samples:
    
```     
> plotHist(mod)
```          
 - _optional_ Produce the library of TCR sequences associated with positve samples:

```     
> getLib(mod) 
```     
5) Predict sample(s) exposure based on model:
```
> pred(mod, "path/to/unknown", "unknown-sample-label", FIELD)
```     
_Note_: If predicting multiple samples, use vectors for paths and labels.
       
       
