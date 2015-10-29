# titanCNAutils

An R package for working TitanCNA Results

To install this package using devtools:

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
