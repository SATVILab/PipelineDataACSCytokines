test_that(".fix_names works", {
  test_list <- list(
    "a", "b",
    "c" = "1",
    "a\\c" = 1, "b/def" = 2, 
    "c" = "3"
  )
  expect_identical(
     names(.fix_names(test_list)),
     c("v1", "v2", "c_1", "a_c", "b_def", "c_2")
  )

})
