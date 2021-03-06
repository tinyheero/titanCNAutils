% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/get_overlap_seg_with_mask_ind.R
\name{get_overlap_seg_with_mask_ind}
\alias{get_overlap_seg_with_mask_ind}
\title{Find Overlapping Segments with a CNV Mask}
\usage{
get_overlap_seg_with_mask_ind(segs.df, mask.df, overlap.prop = 0.25)
}
\arguments{
\item{segs.df}{TitanCNA segments in a data.frame. This can be the
data.frame loaded from the \code{load_titan_seg} function. Or at the very
least the data.frame should contain the columns: Chromosome,
Start_Position, End_Position}

\item{mask.df}{Mask segment in a data.frame. Must contain the columns: chr,
start, end}

\item{overlap.prop}{Proportion that the TitanCNA segment must overlap with
cnv mask segment in order to be considered overlapping}
}
\value{
Vector containing the row number in the segs.df that overlaps with
 any cnv.mask segment
}
\description{
Given a TitanCNA segment data.frame and a mask segment data.frame, this
function will find whether any TitanCNA segment overlaps, with any
significant degree according to the overlap.prop parameter, with a CNV mask
segment. It returns the row number of any segment that does.
}
\examples{
segs.df <- data.frame(chr = c("chr1", "chr2", "chr1", "chr3"),
                      start = 1:4, end = 7:10,
                      stringsAsFactors = FALSE)
mask.df <- data.frame(chr = "chr1", start = 6, end = 8,
                          stringsAsFactors = FALSE)
get_overlap_seg_with_mask_ind(segs.df, mask.df)
}

