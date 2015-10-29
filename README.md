# titanCNAutils

An R package for working [TitanCNA](https://github.com/gavinha/TitanCNA) Results.

As this is not part of Bioconductor yet, you will need to manually install some Bioconductor packages first:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("S4Vectors")
biocLite("IRanges")
```

Once these Bioconductor packages are installed, you can install this package using devtools:

```{r}
devtools::install_github("tinyheero/titanCNAutils")
```

# Overview

To see the full list of exported functions:

```{r}
library("titanCNAutils")
ls("package:titanCNAutils")
```

A quick overview of some of the key functions:

* load_titan_seg: Loads a TitanCNA segment file generated from the createTITANsegmentfiles.pl script
* get_overlap_seg_with_mask_ind: Identify TitanCNA segments that overlap with potential germline CNV segments
* get_gene_titan_seg_overlap: Annotates the TitanCNA segments with gene annotations
