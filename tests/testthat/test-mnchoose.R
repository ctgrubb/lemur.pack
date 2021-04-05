test_that(
  "(l)mnchoose; test examples",
  {
    expect_equal(mnchoose(7, c(3, 2, 2)), 210)
    expect_equal(lmnchoose(7, c(3, 2, 2)), log(210))
  }
)
