# Precompiled vignettes that take too long
# Must manually move image files from lemur.pack/ to lemur.pack/vignettes/ after knit

library(knitr)

system("rm -rf figure/")
knit("vignettes/articles/discreteIntro.Rmd.orig", "vignettes/articles/discreteIntro.Rmd")
system("rm vignettes/articles/figure/*")
system("mv figure/* vignettes/articles/figure/")
system("rm -rf figure/")
