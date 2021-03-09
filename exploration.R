n <- 12
k <- c(3, 4, 5)

test <- factorial(12) / (factorial(3) * factorial(4) * factorial(5))
test2 <- log(test)

sum(log(1:12)) - sum(log(1:3)) - sum(log(1:4)) - sum(log(1:5))

sapply(k, function(x) sum(log(1:x)))

