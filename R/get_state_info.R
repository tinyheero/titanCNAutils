#' Decode the Titan State To Give Copy number and State Name.
#'
#' @param G titan state
#' @param symmetric Boolean flag to indicate whether "similar" states should be collapsed
#' @author Gavin Ha 
#' \url{https://github.com/gavinha/TitanCNA/blob/master/R/utils.R}
decodeLOH <- function(G, symmetric = TRUE) {
  T <- length(G)
  Z <- rep("NA", T)
  CN <- rep(NA, T)

  if (symmetric) {
    DLOH <- G == 1
    NLOH <- G == 2
    ALOH <- G == 4 | G == 6 | G == 9 | G == 12 |
      G == 16 | G == 20
    HET <- G == 3
    GAIN <- G == 5
    ASCNA <- G == 7 | G == 10 | G == 13 | G ==
      17 | G == 21
    BCNA <- G == 8 | G == 15 | G == 24
    UBCNA <- G == 11 | G == 14 | G == 18 | G ==
      19 | G == 22 | G == 23
  } else {
    DLOH <- G == 1 | G == 2
    NLOH <- G == 3 | G == 5
    ALOH <- G == 6 | G == 9 | G == 10 | G == 14 |
      G == 15 | G == 20 | G == 22 | G == 28 |
      G == 29 | G == 36 | G == 37 | G == 45
    HET <- G == 4
    GAIN <- G == 7 | G == 8
    ASCNA <- G == 11 | G == 13 | G == 16 | G ==
      19 | G == 23 | G == 27 | G == 30 | G ==
      35 | G == 38 | G == 44
    BCNA <- G == 12 | G == 25 | G == 41
    UBCNA <- G == 17 | G == 18 | G == 24 | G ==
      26 | G == 31 | G == 32 | G == 33 | G ==
      34 | G == 39 | G == 40 | G == 42 | G ==
      43
  }
  HOMD <- G == 0
  OUT <- G == -1

  Z[HOMD] <- "HOMD"
  Z[DLOH] <- "DLOH"
  Z[NLOH] <- "NLOH"
  Z[ALOH] <- "ALOH"
  Z[HET] <- "HET"
  Z[GAIN] <- "GAIN"
  Z[ASCNA] <- "ASCNA"
  Z[BCNA] <- "BCNA"
  Z[UBCNA] <- "UBCNA"
  Z[OUT] <- "OUT"

  if (symmetric) {
    CN[HOMD] <- 0
    CN[DLOH] <- 1
    CN[G >= 2 & G <= 3] <- 2
    CN[G >= 4 & G <= 5] <- 3
    CN[G >= 6 & G <= 8] <- 4
    CN[G >= 9 & G <= 11] <- 5
    CN[G >= 12 & G <= 15] <- 6
    CN[G >= 16 & G <= 19] <- 7
    CN[G >= 20 & G <= 24] <- 8
  } else {
    CN[HOMD] <- 0
    CN[DLOH] <- 1
    CN[G >= 3 & G <= 5] <- 2
    CN[G >= 6 & G <= 9] <- 3
    CN[G >= 10 & G <= 14] <- 4
    CN[G >= 15 & G <= 20] <- 5
    CN[G >= 21 & G <= 28] <- 6
    CN[G >= 29 & G <= 36] <- 7
    CN[G >= 37 & G <= 45] <- 8
  }

  output <- vector("list", 0)
  output$G <- Z
  output$CN <- CN
  output
}

#' Get Titan State Information
#'
#' Returns information about the Titan states in a list structure
#'
#' @return A list with names corresponding to various details "categories" of
#'   state information
#' @export
get_state_info <- function() {

  my.list <- list(
    "states" = factor(0:24),
    "states.del" = c(0,1),
    "states.del.loh" = c(0:2), # counts deletion and the somatic 2N LOH state 
    "states.neu" = c(3), # counts CN neutral state, but excludes the somatic 2N LOH state
    "states.neu.loh" = c(2,3), # counts CN neutral states including the somatic 2N LOH state
    "states.gain" = c(4:24), # counts all CN gain/amplified states
    "states.somloh" = c(2,4,7,9,12,16,20), # counts all somatic LOH states
    "states.HOMD" = c(0), 
    "states.HETD" = c(1),
    "states.cn.neu" = c(3), # counts CN neutral states that excludes the somatic 2N LOH state
    "states.cn.neu.loh" = c(2), # counts somatic 2N LOH state 
    "states.HLAMP" = c(9:24), # anything higher than 3N
    "states.cols.summary" = c("HOMD" = "#1F78B4", 
                              "HETD" = "#A6CEE3", 
                              "NEU" = "lightgrey",
                              "3N" = "#FB9A99", 
                              "4N" = "#FDAE6B",
                              "5N" = "#F16913",
                              "6N" = "#E31A1C", 
                              "7N" = "#9E1214", 
                              "8N" = "#590A0B",
                              "SOMLOH" = "#33A02C")
    ) 

  titan.states <- as.numeric(as.character(my.list[["states"]]))
  titan.states.decode <- decodeLOH(titan.states)
  titan.states.decode.df <- lapply(titan.states.decode, data.frame)
  titan.states.decode.df <- dplyr::bind_cols(titan.states.decode.df)
  colnames(titan.states.decode.df) <- c("stateName", "copyNum")
  titan.states.decode.df <- dplyr::mutate(
                              titan.states.decode.df, 
                              stateNameModified = paste(copyNum, "N_", stateName, sep = "")
                              )

  my.list[["states.names"]] <- setNames(titan.states.decode.df[, "stateNameModified"], 0:24)

  # this puts the states into overall copy number classes
  my.list[["states.names.summary"]] <- setNames(
                                        c("HOMD", "HETD", "SOMLOH", 
                                          "NEU", "SOMLOH", "3N", "SOMLOH", 
                                          "4N", "4N", "SOMLOH", "5N", "5N", 
                                          "SOMLOH", "6N", "6N", "6N", "SOMLOH", 
                                          "7N", "7N", "7N", "SOMLOH", "8N", 
                                          "8N", "8N", "8N" ), 0:24) 

  my.list[["copy.num"]] <- titan.states.decode.df[, "copyNum"]
  my.list
}
