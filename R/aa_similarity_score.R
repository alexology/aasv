#' Score based on aa similarity
#'
#' @description
#' This function calculates the string similarity, the grantham distance and the phylogenetic
#' conservation of amino acids against a reference database.
#'
#' @param folder_path Path to the folder where the results are stored.
#' @param output_path Path to the results folder. Default to folder path.
#' @param asv_taxonomy A file with the taxonomy assigned to ASV. It must have at least
#' 3 columns: family, genus and species.
#' @param aa_sim_mat Pairwise similarity between triplets calculated with the function
#' \code{aa_similarity_matrix}.
#' @param genetic_code Genetic code as defined in the \code{getGeneticCode} function
#' of the \code{Biostrings} package.
#' @param FUN_sim Function to calculate the similarity among query and reference
#' triplets. Default to \code{min}.
#' @param FUN_gra Function to calculate the grantham distance among query and reference
#' triplets. Default to \code{min}.
#'
#' @details
#' The triplet of an amino acid at a specific position that is not found at the same
#' position in a reference dataset is compared with triplets of the same reference database.
#' Three indices are calculated: string similarity, grantham distance and phylogenetic conservation.
#' String similarity is calculated with the \code{hamming} distance of the
#' \code{stringdist} package. The grantham distance is calculated with the package
#' \code{grantham}, while the phylogenetic conservation is calculated by comparing
#' how much an amino acid is conserved at a given position of the reference database.
#' Continuous results are transformed into classes: conservative, moderately conservative,
#' moderately radical and radical.
#'
#'
#' @return A tibble with the same columns as \code{asv_taxonomy}, joined with
#'   \code{triplet_similarity}, \code{grantham_score}, and
#'   \code{phylogenetic_score} columns from the computed results.
#'
#' @importFrom Biostrings readAAStringSet subseq vmatchPattern
#' @importFrom dplyr bind_rows filter pull rename slice
#' @importFrom grantham as_three_letter grantham_distance
#' @importFrom tibble as_tibble
#' @importFrom utils read.table
#'
#' @export

# https://onlinelibrary.wiley.com/doi/full/10.1002/prot.10458?casa_token=PHnG0u2LJhYAAAAA%3AzW1PvqBLTeY2HftXmsIVlw5kJprmofjk62lsdIUVfwAOtDD-FowaJ2VIauAn9pDrEulaUyV75xVHnZM
# https://academic.oup.com/nar/article/45/3/e13/2972199?login=false


aa_similarity_score <- function(folder_path = NULL,
                                output_path = NULL,
                                asv_taxonomy = NULL,
                                aa_sim_mat = NULL,
                                genetic_code = "5",
                                FUN_sim = min,
                                FUN_gra = min){

  # set the folder path, default to working directory
  if(is.null(folder_path)){
    folder_path <- getwd()
  }

  # set the output path, default to folder path
  if(is.null(output_path)){
    output_path <- folder_path
  }

  # get specie list
  query_taxa <- asv_taxonomy %>%
    dplyr::pull(species) %>%
    unique()

  # set the data.frame to store the results
  res <- data.frame(sequences = character(),
                    tax_lev = character(),
                    error = character(),
                    triplet_similarity = numeric(),
                    grantham_score = numeric(),
                    phylogenetic_score = character())

  for (i in seq_along(query_taxa)) {
    # get the ith species name
    query_taxa_i <- query_taxa[i]

    # get family name
    species_position <- which(asv_taxonomy == query_taxa_i, arr.ind = TRUE)

    # get family name
    family_name <- asv_taxonomy %>%
      dplyr::slice(species_position[1]) %>%
      dplyr::pull(family) %>%
      tibble::as_tibble()

    # iterate through taxonomic levels
    tax_lev <- c("species", "genus", "family")

    # subset asv_taxonomy for further calculations
    asv_taxonomy_i <- asv_taxonomy %>%
      dplyr::filter(species %in% query_taxa_i)

    for (j in seq_along(tax_lev)) {
      # path to alignments
      ali_path <- file.path(output_path,
                            "query_alignments",
                            family_name,
                            tax_lev[j],
                            paste(gsub(" ", "_",query_taxa_i),
                                  ".fasta",
                                  sep = ""))

      # check if file exist, if not skip to the next iteration
      if(! file.exists(ali_path)){
        next()
      } else {
        # get the alignment
        ali_i <- Biostrings::readAAStringSet(ali_path)
      }


      # errors path
      aa_path <- file.path(output_path,
                           "aa_errors",
                           family_name,
                           "species",
                           paste(gsub(" ", "_",query_taxa_i),
                                 ".txt",
                                 sep = ""))


      # check if file exist, if not skip to the next iteration
      if(! file.exists(aa_path)){
        next()
      } else {
        # get the alignment
        aa_i <- utils::read.table(aa_path, tryLogical = FALSE)
      }

      # depending on the taxonomic level, need to load multiple files
      # get the file list
      aa_tbl_list <- list.files(file.path(folder_path,
                                          "aa_tables",
                                          family_name),
                                full.names = TRUE,
                                pattern = "\\.txt$")

      # skip if no files are available
      if(length(aa_tbl_list) < 1){
        next()
      } else{
        # import files
        a_tbl_i <- lapply(aa_tbl_list, read.table) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::filter(seq_id %in% names(ali_i))

        # skip if no sequences are retained
        if(nrow(a_tbl_i) == 0){
          next()
        }

      }




      # # path to aa_tables
      # a_tbl_path <- file.path(folder_path,
      #                         "aa_tables",
      #                          family_name,
      #                          paste(gsub(" ", "_", query_taxa_i),
      #                                ".txt",
      #                                sep = ""))

      # # check if file exist, if not skip to the next iteration
      # if(! file.exists(a_tbl_path)){
      #   next()
      # } else {
      #   # get the alignment, while keeping only sequences in the alignment
      #   a_tbl_i <- utils::read.table(a_tbl_path) %>%
      #     dplyr::filter(seq_id %in% names(ali_i))
      #
      #   if(nrow(a_tbl_i) == 0){
      #     next()
      #   }
      # }


      # iterate through query sequences
      for (z in seq_len(nrow(aa_i))) {

        ### SEQUENCE SIMILARITY ------------------------------------------------
        # get ith query
        query_i <- aa_i %>%
          dplyr::slice(z) %>%
          dplyr::pull("query_sequence")

        # get query number
        query_number <- gsub("query_", "", query_i)

        # get alignments name
        ali_names <- names(ali_i)

        # position of the error in the alignment
        query_sequence <- asv_taxonomy_i %>%
          dplyr::slice(as.numeric(query_number)) %>%
          dplyr::pull(ASV)

        # get wrong aa
        start_pos <- best_translation(query_sequence,
                                                  genetic_code = genetic_code,
                                                  res_type = "position")

        # get the position of the amino acid to check
        to_check <- aa_i %>%
          dplyr::slice(z) %>%
          dplyr::pull(position)

        # get the position of the triplet
        q_i_triplet_start <- start_pos + (to_check - 1) * 3

        # get the triple
        q_i_triplet <- substring(query_sequence,
                                 first = q_i_triplet_start,
                                 last= q_i_triplet_start + 2)

        # get the translation of the triplet
        triplet_translated <- best_translation(q_i_triplet,
                                                            genetic_code = genetic_code)
        # check if the translation is right
        query_aa <-  aa_i %>%
          dplyr::slice(z) %>%
          dplyr::pull(aa)

        # if translation is wrong stop
        if(! identical(query_aa, triplet_translated)){
          stop("Lost in translation.")
        }

        # subset ali
        ali_i_sub <- ali_i[which(query_i %in% names(ali_i))]

        # position of the wrong aminoacid in the alignment
        wrong_position <- ali_i[which(query_i %in% names(ali_i))] %>%
          as.character() %>%
          substring(., first = 1:nchar(.), last = 1:nchar(.))

        # first position of the query sequence i the alignment
        first_position <- which(!grepl("-", wrong_position))[1]

        ### CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # position in the alignment
        ali_position <- first_position + to_check - 1

        ### WARNING!!! TO REVIEW!!!
        # aa_right <- subset(a_tbl_i, position == ali_position) %>%
        #   dplyr::pull(triplet)

        unique_right_names <- unique(a_tbl_i$seq_id)

        aa_right <- c()

        for (k in seq_along(unique_right_names)) {
          aa_temp <- a_tbl_i %>%
            dplyr::filter(seq_id == unique_right_names[k])

          # aligned string for THIS reference (was erroneously re-reading the
          # query's own aligned string on every iteration, so every reference
          # was compared against the query's offset instead of its own)
          ali_i_temp <- ali_i[which(names(ali_i) %in% unique_right_names[k])] %>%
            as.character()

          if(length(ali_i_temp) == 0){
            next()
          }

          # convert the alignment-column index (ali_position) into this
          # reference's own ungapped aa_table position, by counting its
          # non-gap characters up to that column - robust to internal gaps,
          # matching the equivalent idiom used in calculate_asv_distance.R
          real_position <- nchar(gsub("-", "", substring(ali_i_temp, 1, ali_position)))

          triplet <- aa_temp %>%
            dplyr::filter(position == real_position) %>%
            dplyr::pull(triplet)

          aa_right <- c(aa_right, triplet)

        }



        aa_result <- aa_sim_mat %>%
          dplyr::filter(triplets_1 == q_i_triplet) %>%
          dplyr::filter(triplets_2 %in% aa_right) %>%
          dplyr::pull(value) %>%
          FUN_sim()

        ### GRANTHAM DISTANCE --------------------------------------------------

        aa_right_gra <- ali_i[which(names(ali_i) %in% c(unique(a_tbl_i$seq_id)))] %>%
          Biostrings::subseq(., start = real_position, end = real_position) %>%
          grantham::as_three_letter()

        aa_right_gra <- aa_right_gra[!is.na(aa_right_gra)]

        aa_wrong_gra <- query_aa %>%
          grantham::as_three_letter() %>%
          rep(., length = length(aa_right_gra))


        grantham_dist <- grantham::grantham_distance(aa_right_gra, aa_wrong_gra) %>%
          dplyr::pull(d) %>%
          FUN_gra()

        ### CONSERVATION -------------------------------------------------------

        # calculate aa frequency in the reference database
        cons_frequencies <- as.vector(table(aa_right_gra))

        # calculate the frequency of of each amino acid
        cons_frequencies <- max(cons_frequencies)/sum(cons_frequencies)


        # get the results into a data.frame
        res_temp <- data.frame(sequences = query_sequence,
                               tax_lev = tax_lev[j],
                               error = query_i,
                               triplet_similarity = aa_result,
                               grantham_score = grantham_dist,
                               phylogenetic_score = phylogeny_score(cons_frequencies))

        # append the results
        res <- dplyr::bind_rows(res, res_temp)

      }

    }

    message(query_taxa_i, ": done")
  }

  # res %>%
  #   tibble::as_tibble(.) %>%
  #   mutate(triplet_similarity = aaplotype:::triplet_similarity(triplet_similarity)) %>%
  #   mutate(grantham_score = aaplotype:::grantham_score(grantham_score)) %>%
  #   inner_join(taxonomy, ., by = "sequences") %>%
  #   dplyr::rename(sequence = sequences)

  res %>%
    tibble::as_tibble(.) %>%
    # dplyr::mutate(triplet_similarity = triplet_similarity) %>%
    # dplyr::mutate(grantham_score = grantham_score) %>%
    dplyr::inner_join(asv_taxonomy, ., by = "sequences") %>%
    dplyr::rename(sequence = sequences)

}
