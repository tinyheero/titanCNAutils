#' Get Gene Overlaps with TitanCNA segments
#'
#' @param titan.seg.gr GRanges object of the TitanCNA seg
#' @param gene.annot.gr GRanges object of the genes 
#' @return data.frame that contains the overlap of the input gene annotation
#'  with the input segment data.
#' @export
get_gene_titan_seg_overlap <- function(titan.seg.df, gene.annot.df) {

  message("Creating GRanges")
  titan.seg.gr <- 
    GenomicRanges::makeGRangesFromDataFrame(titan.seg.df, 
                                            seqnames.field = "Chromosome", 
                                            start.field = "Start_Position",
                                            end.field = "End_Position", 
                                            keep.extra.columns = TRUE)

  gene.annot.gr <- 
    GenomicRanges::makeGRangesFromDataFrame(gene.annot.df, 
                                            keep.extra.columns = TRUE)

  message("Getting Overlaps")
  overlap.res <- GenomicRanges::findOverlaps(gene.annot.gr, titan.seg.gr)
  gene.annot.gr.hits <- S4Vectors::queryHits(overlap.res)
  titan.seg.gr.hits <- S4Vectors::subjectHits(overlap.res)

  gene.annot.gr.overlap <- gene.annot.gr[gene.annot.gr.hits, ]
  titan.seg.gr.overlap.mcols <- 
    S4Vectors::mcols(titan.seg.gr[titan.seg.gr.hits, ])

  # Build output data.frame
  overlap.df <- BiocGenerics::cbind(S4Vectors::mcols(gene.annot.gr.overlap), 
                                    titan.seg.gr.overlap.mcols)
  overlap.df <- as.data.frame(overlap.df)
  overlap.df <- dplyr::as_data_frame(overlap.df)

  # Ensures that genes with no overlaps also get reported
  out.df <- dplyr::left_join(gene.annot.df, overlap.df)
  out.df
}
