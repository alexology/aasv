# https://stackoverflow.com/questions/12403312/find-the-number-of-spaces-in-a-string
countSpaces <- function(s) {
  sapply(gregexpr(" ", s), function(p) sum(p >= 0))
}

remove_minus <- function(x, forward = TRUE) {
  if (forward) {
    for (i in seq_along(x)) {
      if (x[1] == -1) x <- x[-1] else break
    }
  } else {
    for (i in seq_along(x)) {
      if (x[length(x)] == -1) x <- x[-length(x)] else break
    }
  }
  x
}

replace_minus <- function(x, forward = TRUE) {
  if (forward) {
    for (i in seq_along(x)) {
      if (x[i] == -1) x[i] <- NA else break
    }
  } else {
    for (i in seq_along(x)) {
      j <- length(x) - i + 1
      if (x[j] == -1) x[j] <- NA else break
    }
  }
  x
}

best_translation <- function(x, genetic_code = NULL, res_type = "sequence") {

  i   <- 1
  res <- NA
  max_iter <- nchar(as.character(x))

  while (is.na(res) && i <= max_iter) {
    translation <- x %>%
      Biostrings::DNAStringSet() %>%
      Biostrings::subseq(i) %>%
      Biostrings::translate(genetic.code = Biostrings::getGeneticCode(genetic_code),
                            no.init.codon = TRUE) %>%
      suppressWarnings()

    has_stop <- grepl("*", as.character(translation), fixed = TRUE)

    if (has_stop) {
      i <- i + 1
    } else {
      res <- as.character(translation)
    }
  }

  if (identical(res_type, "sequence")) res else i
}

inner_gaps <- function(x) {
  z <- apply(x, 2, function(col) ifelse(col == "-", -1, 1))

  if (is.vector(z)) {
    z <- t(as.matrix(z))
  }

  logical_store <- logical(nrow(z))

  for (w in seq_len(nrow(z))) {
    z1 <- z[w, ]

    if (z1[1] == -1)          z1 <- remove_minus(z1, forward = TRUE)
    if (z1[length(z1)] == -1) z1 <- remove_minus(z1, forward = FALSE)

    z2 <- z1[-length(z1)] * z1[-1]
    logical_store[w] <- any(z2 == -1)
  }

  logical_store
}

grantham_score <- function(x) {
  grantham_thresh <- c(0, 50, 100, 150, 1000)
  cut(x,
      grantham_thresh,
      labels = c("conservative", "moderately conservative",
                 "moderately radical", "radical"),
      right = TRUE,
      include.lowest = TRUE) %>%
    as.character()
}

triplet_similarity <- function(x) {
  triplet_thresh <- c(0, 0.1, 0.33, 0.67, 10)
  cut(round(x, 2),
      triplet_thresh,
      labels = c("radical", "moderately radical",
                 "moderately conservative", "conservative"),
      right = TRUE,
      include.lowest = TRUE) %>%
    as.character()
}

phylogeny_score <- function(x) {
  ifelse(x == 1, "conservative", "radical")
}

# Produces a codon-level table for a DNAStringSet.
# Returns a data.frame with columns: seq_id, position, triplet, aa.
# Used by align_species_seq to populate the aa_tables output folder.
aa_table <- function(x, genetic_code = "5") {
  gc <- Biostrings::getGeneticCode(genetic_code)

  result <- lapply(seq_along(x), function(i) {
    seq_name <- names(x)[i]
    dna_str  <- as.character(x[[i]])
    start    <- best_translation(as.character(x[i]),
                                 genetic_code = genetic_code,
                                 res_type     = "position")

    trimmed  <- substring(dna_str, start)
    n_codons <- nchar(trimmed) %/% 3L
    if (n_codons == 0L) return(NULL)

    starts  <- seq(1L, n_codons * 3L - 2L, by = 3L)
    ends    <- seq(3L, n_codons * 3L,       by = 3L)
    triplets <- substring(trimmed, starts, ends)
    aas      <- vapply(triplets,
                       function(t) if (nchar(t) == 3L) as.character(gc[t]) else NA_character_,
                       character(1L))

    data.frame(seq_id   = seq_name,
               position = seq_along(triplets),
               triplet  = triplets,
               aa       = aas,
               stringsAsFactors = FALSE)
  })

  do.call(rbind, Filter(Negate(is.null), result))
}
