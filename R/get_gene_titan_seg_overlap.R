#' Get Gene Overlaps with TitanCNA segments
#'
#' @param segs.df data.frame of the TitanCNA segs file loaded from 
#'  \code{load_titan_seg}
#' @param gene.annot.df data.frame of the gene annotations
#' @param overlap.prop Proportion that the TitanCNA segment must overlap with
#'   a gene in order to be considered overlapping. The smaller the value the 
#'   more focal the copy number alteration you want to associate with gene
#' @return data.frame that contains the overlap of the input gene annotation
#'  with the input segment data.
#' @export
get_gene_titan_seg_overlap <- function(segs.df, gene.annot.df, overlap.prop = 0.5) {

  if (overlap.prop > 1) {
    stop("overlap.prop cannot be greater than 1")
  }

  message("Creating GRanges")
  titan.seg.gr <- 
    GenomicRanges::makeGRangesFromDataFrame(segs.df, 
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

  message(paste("Checking if segment overlaps are >=", overlap.prop, 
                "of genes"))
  # See if the overlap is at least `overlap.prop` of the segment length
  gene.annot.overlap.prop <- 
    IRanges::width(IRanges::pintersect(gene.annot.gr[gene.annot.gr.hits], 
                                       titan.seg.gr[titan.seg.gr.hits])) / 
    IRanges::width(gene.annot.gr[gene.annot.gr.hits])

  gene.annot.overlap.ind <- which(gene.annot.overlap.prop >= overlap.prop)
  overlap.df <- overlap.df[gene.annot.overlap.ind, ]

  # Ensures that genes with no overlaps also get reported
  out.df <- dplyr::left_join(gene.annot.df, overlap.df)
  out.df
}
