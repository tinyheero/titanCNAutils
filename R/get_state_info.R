#' Decode the Titan State To Give Copy number and State Name.
#'
#' @param G titan state
#' @param symmetric Boolean flag to indicate whether "similar" states should be collapsed
#' @author Gavin Ha 
#' \url{https://github.com/gavinha/TitanCNA/blob/master/R/utils.R}
#' @export
decode_LOH <- function(G, symmetric = TRUE) {
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
#' @return A named list containing various details "categories" of state 
#'   information
#' @export
get_state_info <- function() {

  my.list <- list(
    "states" = factor(0:24),
    "states.del" = c(0,1),

    # counts deletion and the somatic 2N LOH state 
    "states.del.loh" = c(0:2), 

    # counts CN neutral state, but excludes the somatic 2N LOH state
    "states.neu" = c(3), 

    # counts CN neutral states including the somatic 2N LOH state
    "states.neu.loh" = c(2,3),

    # counts all CN gain/amplified states
    "states.gain" = c(4:24),

    # counts all somatic LOH states
    "states.somloh" = c(2,4,7,9,12,16,20),
    "states.HOMD" = c(0), 
    "states.HETD" = c(1),

    # counts CN neutral states that excludes the somatic 2N LOH state
    "states.cn.neu" = c(3), 

    # counts somatic 2N LOH state 
    "states.cn.neu.loh" = c(2), 

    # anything higher than 3N
    "states.HLAMP" = c(9:24),
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

  # Used to assess the "importance" of TitanCNA states
  titan.states.order <- 
    c(-1, # outlier state (OUT)
      3, # 2N copy neutral heterozygous (AB; HET; 2 CN)
      2, # 2N copy neutral LOH (AA/AB; NLOH; 2 CN)
      5, # 3N Allele specific gain (AAB/ABB; GAIN; 3 CN)
      4, # 3N Gain LOH (AAA/BBB; ALOH; 3 CN)
      8, # 4N balanced gain (AABB; BCNA; 4 CN)
      7, # 4N allele specific gain (AAAB/ABBB; ASCNA; 4 CN)
      6, # 4N amplified LOH (AAAA/BBBB; ALOH; 4 CN)
      11, # 5N unbalanced amplification (AAABB/AABBB; UBCNA; 5 CN)
      10, # 5N allele specific amplification (AAAAB/ABBBB; ASCNA; 5 CN)
      9, # 5N amplified LOH (AAAAA/BBBBB; ALOH; 5 CN)
      15, # 6N balanced amplification (AAABBB; BCNA; 6 CN)
      14, # 6N unbalanced amplification (AAAABB/AABBBB; UBCNA; 6 CN)
      13, # 6N allele specific amplification (AAAAAB/ABBBBB; ASCNA; 6 CN)
      12, # 6N amplified LOH (AAAAAA/BBBBBB; ALOH; 6 CN)
      18, # 7N unbalanced amplification (AAAAABB/AABBBBB; UBCNA; 7 CN)
      19, # 7N unbalanced amplification (AAABBBB/AAAABBB; UBCNA; 7 CN)
      17, # 7N allele specific amplifcation (AAAAAAB/ABBBBBB; ASCNA; 7 CN)
      16, # 7N amplified LOH (AAAAAAA/BBBBBBB; ALOH; 7 CN)
      24, # 8N balanced amplification (AAAABBBB; BCNA; 8 CN)
      22, # 8N unbalanced amplification (AAAAAABB/AABBBBBB; UBCNA; 8 CN)
      23, # 8N unbalanced amplification (AAAAABBB/AAABBBBB; UBCNA; 8 CN)
      21, # 8N allele specific amplifcation (AAAAAAAB/ABBBBBBB; ASCNA; 8 CN)
      20, # 8N amplified LOH (AAAAAAAA/BBBBBBBBB; ALOH; 8 CN)
      1, # 1N heterozygous deletion (A/B; DLOH; 1 CN)
      0 # 0N homozygous deletion (HOMD; 0 CN)
		)

  my.list[["titan.states.order"]] <- titan.states.order

  titan.states <- as.numeric(as.character(my.list[["states"]]))
  titan.states.decode <- decode_LOH(titan.states)
  titan.states.decode.df <- lapply(titan.states.decode, data.frame)
  titan.states.decode.df <- dplyr::bind_cols(titan.states.decode.df)
  colnames(titan.states.decode.df) <- c("stateName", "copyNum")
  titan.states.decode.df <- 
    dplyr::mutate(titan.states.decode.df, 
                  stateNameModified = paste(copyNum, "N_", stateName, sep = ""))

  my.list[["states.names"]] <- 
    setNames(titan.states.decode.df[["stateNameModified"]], 0:24)

  # Puts the States into Overall Copy Number Classes
  my.list[["states.names.summary"]] <- 
    setNames(c("HOMD", "HETD", "SOMLOH", "NEU", "SOMLOH", "3N", "SOMLOH", 
               "4N", "4N", "SOMLOH", "5N", "5N", "SOMLOH", "6N", "6N", "6N", 
               "SOMLOH", "7N", "7N", "7N", "SOMLOH", "8N", "8N", "8N", "8N"), 
             0:24) 

  my.list[["copy.num"]] <- 
    setNames(titan.states.decode.df[["copyNum"]], 0:24)
  my.list
}
