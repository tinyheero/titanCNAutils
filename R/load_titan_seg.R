#' Load the TitanCNA Seg File
#'
#' @param titan.seg.file Path to TitanCNA seg file. This should be the file 
#'  generated from the createTITANsegmentfiles.pl script provided in the 
#'  TitanCNA download
#' @return data.frame of the TitanCNA seg file
#' @export
load_titan_seg <- function(titan.seg.file) {
  titan.seg.df <- 
    readr::read_tsv(titan.seg.file, 
                    col_types = list("Chromosome" = readr::col_character()))
  titan.seg.df
}

