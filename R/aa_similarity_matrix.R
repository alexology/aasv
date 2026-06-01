#' Triplet similarity
#'
#' @description
#' This function calculates the pairwise similarity between triplets.
#'
#' @param genetic_code Genetic code as defined in the \code{getGeneticCode} function
#' of the \code{Biostrings} package.
#'
#' @details
#' There is no need to run this function each time. The results can be stored and
#' uploaded when needed.
#'
#' @return A data frame with three columns: \code{triplets_1}, \code{triplets_2},
#'   and \code{value} (pairwise Hamming similarity in [0, 1]).
#'
#' @export
#'
#' @importFrom Biostrings getGeneticCode
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom dplyr distinct_all inner_join mutate pull rename select
#' @importFrom tibble as_tibble rowid_to_column rownames_to_column
#' @importFrom stringdist stringsimmatrix
#'
#' @examples
#' \donttest{
#' aa_sim <- aa_similarity_matrix(genetic_code = "5")
#' }

aa_similarity_matrix <- function(genetic_code = "5"){

  # get the genetic code from Biostrings
  g_code_aa <- Biostrings::getGeneticCode(genetic_code)

  # get names from the vector
  # Basically names are nucleotide triplets whil the lements of the vectors are
  # aminoacids
  g_code_aa_names <- names(g_code_aa)

  # revert, now we have triplets as elements and aminoacids as names
  names(g_code_aa_names) <- g_code_aa

  # create a data.frame with an ID and the triplets
  triplets_df <- data.frame(ID = 1:length(g_code_aa_names),
                            sequence = g_code_aa_names)

  # expand grid to have all the combinations of triplets
  triplets_expanded_df <- tidyr::expand_grid(triplets_1 = triplets_df[, "sequence"],
                                             triplets_2 = triplets_df[, "sequence"]) %>%
    dplyr::distinct_all()

  # get a data.frame for further inner_
  z_aa2 <- triplets_df %>%
    dplyr::select(-1) %>%
    tibble::rowid_to_column("ID")

  # get the first set of triplet as a vector
  triplets_1 <- triplets_expanded_df %>%
    dplyr::pull("triplets_1")

  # get the second set of triplet as a vector
  triplets_2 <- triplets_expanded_df %>%
    dplyr::pull("triplets_2")

  # make triplets_2 a data.frame for a later inner_join
  triplets_2_df <- triplets_2 %>%
    as.data.frame() %>%
    dplyr::rename(triplets_2 = 1) %>%
    tibble::rowid_to_column("ID")

  # calculate similarity with the function stringsimmatrix of the stringdist package
  # and manipulate the results to have a data.frame for further calculations
  aa_sim_mat <- stringdist::stringsimmatrix(triplets_1, triplets_2, method = "hamming") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("triplets_1") %>%
    dplyr::mutate(triplets_1 = triplets_expanded_df$triplets_1) %>%
    tidyr::pivot_longer(-triplets_1) %>%
    dplyr::mutate(name = as.numeric(gsub("V", "", name))) %>%
    dplyr::rename(ID = name) %>%
    dplyr::inner_join(., triplets_2_df, by = "ID") %>%
    dplyr::select(-ID) %>%
    dplyr::select(triplets_1, triplets_2, value) %>%
    dplyr::distinct_all()


  # return the results
  aa_sim_mat
}
