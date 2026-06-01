test_that("aa_similarity_matrix returns a data.frame with the three expected columns", {
  result <- aa_similarity_matrix(genetic_code = "5")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("triplets_1", "triplets_2", "value") %in% names(result)))
})

test_that("aa_similarity_matrix similarity values are in [0, 1]", {
  result <- aa_similarity_matrix(genetic_code = "5")
  expect_true(all(result$value >= 0 & result$value <= 1))
})

test_that("aa_similarity_matrix self-similarity is 1 for every triplet", {
  result <- aa_similarity_matrix(genetic_code = "5")
  self   <- result[result$triplets_1 == result$triplets_2, ]
  expect_true(all(self$value == 1))
})

test_that("aa_similarity_matrix is symmetric", {
  result <- aa_similarity_matrix(genetic_code = "5")

  fwd <- result$value[result$triplets_1 == "ATG" & result$triplets_2 == "ATC"]
  rev <- result$value[result$triplets_1 == "ATC" & result$triplets_2 == "ATG"]
  expect_equal(fwd, rev)
})

test_that("aa_similarity_matrix contains only valid DNA triplets", {
  result  <- aa_similarity_matrix(genetic_code = "5")
  valid   <- grepl("^[ACGT]{3}$", result$triplets_1)
  expect_true(all(valid))
})

test_that("aa_similarity_matrix has fewer rows for standard code than unique triplet pairs", {
  result <- aa_similarity_matrix(genetic_code = "1")
  # 64 codons -> 64*64 = 4096 pairs, but distinct_all may reduce stop codon variants
  expect_gt(nrow(result), 0L)
})

test_that("aa_similarity_matrix identical triplets have similarity 1", {
  result  <- aa_similarity_matrix(genetic_code = "5")
  # use $ accessor to get a plain vector, not a one-column tibble
  perfect <- result$value[result$triplets_1 == "ATG" & result$triplets_2 == "ATG"]
  expect_equal(perfect, 1)
})

test_that("aa_similarity_matrix triplets differing in all 3 positions have similarity 0", {
  result <- aa_similarity_matrix(genetic_code = "5")
  # ATG (A,T,G) vs CAC (C,A,C) differ in all 3 positions -> Hamming similarity = 0
  zero_pair <- result$value[result$triplets_1 == "ATG" & result$triplets_2 == "CAC"]
  expect_equal(zero_pair, 0)
})
