# larimaR
L1000/CMAP Data Exploration and Analysis Tools

LarimaR is designed to help explore and use the L1000 data provided by clue.io

# Demo

## Download data

The first step is to download the CMAP data in gctx format along with the metadata files needed for the extraction.

The data can be located here:

https://clue.io/data/CMap2020#LINCS2020

Download the following files:

* level5_beta_trt_cp_n720216x12328.gctx - expression data for compounds
* cellinfo_beta.txt - cell information
* compoundinfo_beta.txt - compound information
* siginfo_beta.txt - experimental information


## Search for cells/pertunagens with shiny apps

LarimaR comes with 2 shiny based filter apps based on cell and perturbagen data

runCmapCellFilter(cellFile) and runCmapPertFilter(pertFile)

These functions use the files downloaded above, and can be useful to choose a pertubagen and cell of interest.

```{}
dataDir="/home/bj/Research/data"
cellFile="cellinfo_beta.txt"
runCmapCellFilter(cellFile)
```



![cellfilter](https://github.com/user-attachments/assets/7717015e-494e-4fcb-9a64-5f9ce1fd0ef7)

Using this, we can search for lung related cells and export a csv

Similarly, we can use the pertunagen seaech tool to find mtor inhibitors and export a csv of perturbagens.


## Export data

The main function in larimaR is getDosagePipeline() which extracts data from a LINCS compatible dataset.

This function takes in a number of parameters:

* pertName: the string form of a perturbagen
* sigInfo: the data object containing the experimental metadata from a LICS dataset
* gctxFileLocation: the location on disk of the gctx formated data file
* cell: the name of the cell to extract data from
* curphase: the phase of the data used in the extraction. This is needed to resolve changes across phases.

Then we can set our file locations:

```{}
dataDir="/home/bj/Research/data"
sigInfo3=read.csv(paste0(dataDir,"/siginfo_beta.txt"),sep="\t")
compounds=read.csv(paste0(dataDir,"/compoundinfo_beta.txt"),sep="\t")
cells=read.csv(paste0(dataDir,"/cellinfo_beta.txt"),sep="\t")
gctxFile3=paste0(dataDir,"/level5_beta_trt_cp_n720216x12328.gctx")
```

LarimaR has tools to look up internal identfiers using common names:

```{
library(DT)
brds=getBRDS("sirolimus",sigInfo3)
print(brds)
```

 [1] "BRD-K84937637" "BRD-A79768653" "BRD-A50287119" "BRD-A23770159"
 [5] "BRD-K89626439" "BRD-K99369265"

 
