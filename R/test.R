# test_v <- cbind(c(1, -exp(1)), c(-exp(1), 10))
# test_h <- solve(test_v)

# bivariate normal
# mvtnorm::dmvnorm(x, mu, test_v, log=TRUE)
# test_f <- function(x, mu = c(2, 1), v = test_v, h = test_h) {-log(2 * pi) - log(det(v)) / 2 + (-t(x - mu) %*% h %*% (x - mu) / 2)}
# test_g <- function(x, mu = c(2, 1), v = test_v, h = test_h) {as.numeric(-t(x - mu) %*% h)}

# foo <- NUTS::NUTS(x, test_f, test_g, 1024, NULL,         0)
# bar <- NUTS::NUTS(x, test_f, test_g, 1024, diag(test_h), 0)
# baz <- NUTS::NUTS(x, test_f, test_g, 1024, test_h, 0)
# baz <- NUTS::NUTS(x, test_f, test_g, 1024, test_h, 0, eps = 1, find=FALSE)
