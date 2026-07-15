# Factory: returns a performer closure with shared exponential-backoff state.
# Each failure increments the level; each success decays it by one (not a full
# reset), so the backoff stays elevated while the server is struggling and
# gradually drops back as requests start succeeding again.
# Permanent client errors (4xx, other than 408/429 which are transient) are
# not retried: they signal a bad request (e.g. an invalid taxon query) that
# will never succeed, so they are re-raised immediately instead of retrying
# forever.
# Interrupt with Ctrl+C to abort.
make_retry_performer <- function(initial_wait = 60, max_wait = 7200) {
  state <- new.env(parent = emptyenv())
  state$level <- 0L

  function(req, path = NULL) {
    repeat {
      state$level <- state$level + 1L
      if (!is.null(path) && file.exists(path)) unlink(path)
      result <- tryCatch(
        httr2::req_perform(req, path = path),
        error = function(e) {
          status <- e$status
          if (!is.null(status) && status >= 400L && status < 500L &&
              !(status %in% c(408L, 429L))) {
            stop(e)
          }
          NULL
        }
      )
      if (!is.null(result)) {
        state$level <- max(0L, state$level - 1L)
        return(result)
      }
      pause    <- min(initial_wait * 2^(state$level - 1L), max_wait)
      wait_msg <- if (pause >= 60) paste0(round(pause / 60, 1), " min") else paste0(pause, "s")
      message("Request failed (attempt ", state$level, "), retrying in ", wait_msg, "…")
      Sys.sleep(pause)
    }
  }
}

req_perform_retry <- function(req, initial_wait = 60, max_wait = 7200) {
  make_retry_performer(initial_wait, max_wait)(req)
}

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
  # Only the 3 true reading frames are biologically valid starting points.
  # Capping the search here (rather than scanning nucleotide-by-nucleotide
  # to the end of the sequence) means a sequence with a genuine internal
  # stop codon in all 3 frames returns NA - a detectable failure - instead
  # of silently succeeding on some short stop-free tail fragment found by
  # drifting past the stop one nucleotide at a time.
  max_iter <- min(3L, nchar(as.character(x)))

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

  if (identical(res_type, "sequence")) res else if (is.na(res)) NA_integer_ else i
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

causes_internal_gaps <- function(mat) {
  nc      <- ncol(mat)
  non_gap <- mat != "-"
  seq_starts <- apply(non_gap, 1L, function(r) {
    w <- which(r)
    if (length(w)) w[1L] else nc + 1L
  })
  seq_ends <- apply(non_gap, 1L, function(r) {
    w <- which(r)
    if (length(w)) w[length(w)] else 0L
  })
  cols     <- seq_len(nc)
  spanning <- outer(seq_starts, cols, "<=") &
              outer(seq_ends,   cols, ">=")
  n_spanning        <- colSums(spanning)
  n_internal_gap    <- colSums(spanning & !non_gap)
  internal_gap_frac <- ifelse(n_spanning > 1L, n_internal_gap / n_spanning, 0)
  apply(non_gap, 1L, function(r) any(r & (internal_gap_frac > 0.5)))
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
    if (is.na(start)) return(NULL)

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
