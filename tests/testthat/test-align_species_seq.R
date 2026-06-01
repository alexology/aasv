test_that("align_species_seq errors when taxonomy is NULL", {
  expect_error(
    align_species_seq(taxonomy = NULL),
    "`taxonomy` must be a data.frame"
  )
})

test_that("align_species_seq errors when required columns are missing", {
  expect_error(
    align_species_seq(taxonomy = data.frame(family = "x")),
    "`taxonomy` must be a data.frame"
  )
  expect_error(
    align_species_seq(taxonomy = data.frame(family = "x", genus = "y")),
    "`taxonomy` must be a data.frame"
  )
})

test_that("align_species_seq errors when taxonomy column names are wrong", {
  expect_error(
    align_species_seq(taxonomy = data.frame(Family = "x", Genus = "y", Species = "z")),
    "`taxonomy` must be a data.frame"
  )
})
