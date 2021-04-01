library(testit)

assert(
  "lmnchoose(7, c(3, 2, 2)) should be (very close to) log(210)",
  (abs(lmnchoose(7, c(3, 2, 2)) - log(210)) <= sqrt(.Machine$double.eps))
)
