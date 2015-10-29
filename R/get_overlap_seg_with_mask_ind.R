#' Find Overlapping Segments with a CNV Mask
#'
#' Given a TitanCNA segment data.frame and a CNV mask segment, this function
#' will find whether any TitanCNA segment overlaps, with any significant degree
#' according to the overlap.prop parameter, with a CNV mask segment. It returns
#' the row number of any segment that does.
#'
#' @param segs.df TitanCNA segments in a data.frame. This should be the 
#'  data.frame loaded from the load_titan_seg function
#' @param cnv.mask.df CNV mask in a data.frame. Must contain the columns: chr,
#'   start, end
#' @param overlap.prop Proportion that the TitanCNA segment must overlap with
#'   cnv mask segment in order to be considered overlapping
#' @return Vector containing the row number in the segs.df that overlaps with 
#'  any cnv.mask segment
#' @export
#' @examples
#' segs.df <- data.frame(chr = c("chr1", "chr2", "chr1", "chr3"),
#'                       start = 1:4, end = 7:10, 
#'                       stringsAsFactors = FALSE)
#' cnv.mask.df <- data.frame(chr = "chr1", start = 6, end = 8,
#'                           stringsAsFactors = FALSE)
#' get_overlap_seg_with_mask_ind(segs.df, cnv.mask.df)
get_overlap_seg_with_mask_ind <- function(segs.df, cnv.mask.df, 
                                          overlap.prop = 0.25) {

  if (overlap.prop > 1) {
    stop("overlap.prop cannot be greater than 1")
  }

  # Build GRanges Objects
  segs.gr <- 
    GenomicRanges::makeGRangesFromDataFrame(segs.df, 
                                            seqnames.field = "Chromosome", 
                                            start.field = "Start_Position",
                                            end.field = "End_Position")

  cnv.mask.gr <- 
    GenomicRanges::makeGRangesFromDataFrame(cnv.mask.df,
                                            seqnames.field = "chr", 
                                            start.field = "start",
                                            end.field = "end")

  # Find Overlaps
  # The nrows will be larger than the input data since a seg can overlap with 
  # multiple cnv segs. But the values will be indices to the original input 
  # data. This will be used later down the function.
  overlap.res <- GenomicRanges::findOverlaps(segs.gr, cnv.mask.gr)
  gr.hit <- S4Vectors::queryHits(overlap.res)
  subject.hit <- S4Vectors::subjectHits(overlap.res)

  # See if the overlap is at least `overlap.prop` of the segment length
  segs.cnv.overlap.prop <- 
    IRanges::width(IRanges::pintersect(segs.gr[gr.hit], 
                                       cnv.mask.gr[subject.hit])) / 
    IRanges::width(segs.gr[gr.hit])

  segs.overlap.cnv.mask.ind <- which(segs.cnv.overlap.prop >= overlap.prop)

  # Map the indices to the original input data
  unique(gr.hit[segs.overlap.cnv.mask.ind])
}
