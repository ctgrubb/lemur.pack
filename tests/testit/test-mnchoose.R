library(testit)

assert(
  "mnchoose(7, c(3, 2, 2)) should be 210",
  abs(mnchoose(7, c(3, 2, 2)) - 210) < sqrt(.Machine$double.eps)
)
