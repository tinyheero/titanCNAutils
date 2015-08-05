#' Melts the Titan Gene Call Matrix into Tidy Format
#' 
#' This function melts the Titan gene call matrix into a tidy format. It also
#' adds an additional column called stateSummary. This column represents the
#' the states summarized into broad states as specified in the titan.state.info
#' variable
#'
#' @param in.gene.call.matrix The gene call matrix
#' @param titan.state.info Output from the get_state_info() function
#' @return A 4 column data.frame of the Titan gene call data.
#' @export
melt_gene_call_mat <- function(in.gene.call.mat, titan.state.info) {

  out.df <- reshape2::melt(in.gene.call.mat, value.name = "state")
  out.df[, "stateSummary"] <- plyr::revalue(as.character(out.df$state),
                                titan.state.info$states.names.summary)
  out.df
}
