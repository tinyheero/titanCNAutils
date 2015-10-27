#' Find Overlapping Segments with a CNV Mask
#'
#' Given a TitanCNA segment data.frame and a CNV mask segment, this function
#' will find whether any TitanCNA segment overlaps, with any significant degree
#' according to the overlap.prop parameter, with a CNV mask segment. It returns
#' the row number of any segment that does.
#'
#' @param segs.df Titan segments data.frame
#' @param cnv.mask.df CNV mask in a data.frame. Must contain the columns: chr,
#'   start, end
#' @param overlap.prop Proportion that the TitanCNA segment must overlap with
#'   cnv mask segment in order to be considered overlapping
#' @return Vector containing the row number in the segs.df that overlaps with 
#'  any cnv.mask segment
get_overlap_seg_with_mask_ind <- function(segs.df, cnv.mask.df, 
                                          overlap.prop = 0.25) {

  if (overlap.prop > 1) {
    stop("overlap.prop cannot be greater than 1")
  }

  # Build GRanges Objects
  segs.gr <- 
    GenomicRanges::GRanges(seqnames = segs.df$Chromosome,
                           ranges = IRanges::IRanges(
                                      start = segs.df$Start_Position,
                                      end = segs.df$End_Position))
  cnv.mask.gr <- cnv.mask.df %$%
    GenomicRanges::GRanges(seqnames = cnv.mask.df$chr,
                           ranges = IRanges::IRanges(start = cnv.mask.df$start, 
                                                     end = cnv.mask.df$end))

  # Find Overlaps
  overlap.res <- GenomicRanges::findOverlaps(segs.gr, cnv.mask.gr)
  gr.hit <- S4Vectors::queryHits(overlap.res)
  subject.hit <- S4Vectors::subjectHits(overlap.res)

  # See if the overlap is at least `overlap.prop` of the segment length
  segs.cnv.overlap.prop <- 
    IRanges::width(IRanges::pintersect(segs.gr[gr.hit], 
                                       cnv.mask.gr[subject.hit])) / 
    IRanges::width(segs.gr[gr.hit])

  segs.overlap.cnv.mask.ind <- which(segs.cnv.overlap.prop > overlap.prop)
  segs.overlap.cnv.mask.ind
}
