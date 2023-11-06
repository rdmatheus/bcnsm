# The Box-Cox Laplace distribution ------------------------------------------------------------------

#' @name bcla
#'
#' @title The Box-Cox Laplace Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the Box-Cox Laplace distribution with parameters \code{mu},
#'     \code{sigma}, and \code{lambda}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lambda vector of real-valued skewness parameters. If \code{lambda = 0}, the Box-Cox
#'     Laplace distribution reduces to the log-Laplace distribution with parameters
#'     \code{mu} and \code{sigma} (see \code{\link{lla}}).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... further arguments.
#'
#' @return \code{dbcla} returns the density, \code{pbcla} gives the distribution function,
#'     \code{qbcla} gives the quantile function, and \code{rbcla} generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}, with an warning.
#'
#'     The length of the result is determined by \code{n} for \code{rbcla}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'     applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' mu <- 10
#' sigma <- 0.5
#' lambda <- 4
#'
#' # Sample generation
#' x <- rbcla(10000, mu, sigma, lambda)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Box-Cox Laplace Distribution", col = "white")
#' curve(dbcla(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2)
#' legend("topleft", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Box-Cox Laplace Distribution", ylab = "Distribution function")
#' curve(pbcla(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2)
#' legend("topleft", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Box-Cox Laplace Distribution"
#' )
#' curve(qbcla(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL


# Density
#' @rdname bcla
#' @export
dbcla <- function(x, mu, sigma, lambda, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0)] <- NaN

  # Positive density index
  id1 <- which(x > 0 & lambda != 0 & !is.nan(pmf))
  id2 <- which(x > 0 & lambda == 0 & !is.nan(pmf))

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  z[id1] <- ((x[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(x[id2] / mu[id2]) / sigma[id2]

  pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + rLA(z[id1]^2, log = TRUE) -
    RLA(1 / (sigma[id1] * abs(lambda[id1])), log.p = TRUE) -
    lambda[id1] * log(mu[id1]) - log(sigma[id1])

  pmf[id2] <- rLA(z[id2]^2, log = TRUE) - log(sigma[id2] * x[id2])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname bcla
#' @export
pbcla <- function(q, mu, sigma, lambda, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda != 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0)

  z[id1] <- ((q[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(q[id2] / mu[id2]) / sigma[id2]

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda <= 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id1] <- RLA(z[id1]) / RLA(1 / (sigma[id1] * abs(lambda[id1])))
  cdf[id2] <- (RLA(z[id2]) - RLA(-1 / (sigma[id2] * lambda[id2]))) /
    RLA(1 / (sigma[id2] * lambda[id2]))

  cdf[which(q <= 0 & mu > 0 & sigma > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname bcla
#' @export
qbcla <- function(p, mu, sigma, lambda, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(lambda)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda <= 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda > 0)

  zp[id1] <- qLA(p[id1] * RLA(1 / (sigma[id1] * abs(lambda[id1]))))
  zp[id2] <- qLA(1 - (1 - p[id2]) * RLA(1 / (sigma[id2] * abs(lambda[id2]))))

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0)
  id3 <- which(p == 0 & mu > 0 & sigma > 0)
  id4 <- which(p == 1 & mu > 0 & sigma > 0)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * zp[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * zp[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname bcla
#' @export
rbcla <- function(n, mu, sigma, lambda) {
  u <- stats::runif(n)
  qbcla(u, mu, sigma, lambda)
}

# BCS class
bcla <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "bcla"

  # Name
  out$name <- "Box-Cox Laplace"

  # Number of parameters
  out$npar <- 3

  # Extra parameter
  out$extrap <- FALSE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCPE(mu.link = "log"),
                                                      trace = FALSE, tau.fix = TRUE, tau.start = 1), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")),
        exp(stats::coef(gamlss_fit, "sigma")),
        stats::coef(gamlss_fit, "nu"))

    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)
      c(stats::median(x), asinh(CV / 1.5) * gamlss.dist::qPE(0.75, 0, 1, 1), 0L)
    }

  }


  structure(out, class = "bcs")
}

# The log-Laplace distribution ----------------------------------------------------------------------

#' @name lla
#'
#' @title The Log-Laplace Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the log-Laplace distribution with parameters \code{mu} and
#'     \code{sigma}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... further arguments.
#'
#' @details A random variable X has a log-Laplace distribution with parameter \code{mu} and
#'     \code{sigma} if log(X) follows a Laplace distribution with location parameter \code{log(mu)}
#'     and dispersion parameter \code{sigma}. It can be showed that \code{mu} is the median of X.
#'
#' @return \code{dlca} returns the density, \code{plca} gives the distribution
#'     function, \code{qlca} gives the quantile function, and \code{rlca}
#'     generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}.
#'
#'     The length of the result is determined by \code{n} for \code{rlca}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' mu <- 5
#' sigma <- 0.5
#'
#' # Sample generation
#' x <- rlca(500, mu, sigma)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Log-Laplace Distribution", col = "white")
#' curve(dlca(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
#' legend("topright", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Log-Laplace Distribution", ylab = "Distribution function")
#' curve(plca(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
#' legend("bottomright", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Log-Laplace Distribution"
#' )
#' curve(qlca(x, mu, sigma), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL


# Density
#' @rdname lla
#' @export
dlla <- function(x, mu, sigma, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0)] <- NaN

  # Positive density index
  id <- which(x > 0 & !is.nan(pmf))

  # Transformed variables
  z <- rep(NaN, length.out = maxl)

  z[id] <- log(x[id] / mu[id]) / sigma[id]

  pmf[id] <- rLA(z[id]^2, log = TRUE) - log(sigma[id] * x[id])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname lla
#' @export
plla <- function(q, mu, sigma, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  # Transformed variables
  z <- rep(NaN, length.out = maxl)

  id <- which(q > 0 & mu > 0 & sigma > 0)

  z[id] <- log(q[id] / mu[id]) / sigma[id]

  cdf <- rep(NaN, length.out = maxl)
  cdf[id] <- RLA(z[id])

  cdf[which(q <= 0 & mu > 0 & sigma > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname lla
#' @export
qlla <- function(p, mu, sigma, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id <- which(p > 0 & p < 1 & mu > 0 & sigma > 0)
  zp[id] <- qLA(p[id])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0)
  id2 <- which(p == 0 & mu > 0 & sigma > 0)
  id3 <- which(p == 1 & mu > 0 & sigma > 0)

  qtf[id1] <- exp(log(mu[id1]) + sigma[id1] * zp[id1])
  qtf[id2] <- 0
  qtf[id3] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname lla
#' @export
rlla <- function(n, mu, sigma) {
  u <- stats::runif(n)
  exp(log(mu) + sigma * qLA(u))
}

# BCS class
lla <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "lla"

  # Name
  out$name <- "Log-Laplace"

  # Number of parameters
  out$npar <- 2

  # Extra parameter
  out$extrap <- FALSE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCPE(mu.link = "log"),
                                                      trace = FALSE,
                                                      nu.fix = TRUE, nu.fix = 0,
                                                      tau.fix = TRUE, tau.start = 1), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")),
        exp(stats::coef(gamlss_fit, "sigma")))

    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)
      c(stats::median(x), asinh(CV / 1.5) * gamlss.dist::qPE(0.75, 0, 1, 1))
    }

  }

  structure(out, class = "bcs")
}


## Standard Laplace distribution -------------------------------------------------------------------

rLA <- function(u, log = FALSE) {
  pmf <- rep(NaN, length.out = length(u))

  # Positive density index
  id <- which(u >= 0)

  pmf[id] <- -sqrt(u[id]) - log(2)

  if (!log) pmf <- exp(pmf)
  pmf
}

RLA <- function(q, log.p = FALSE) {
  cdf <- rep(0, length.out = length(q))

  id1 <- which(q <= 0)
  id2 <- which(q > 0)

  cdf[id1] <- log(0.5) + q[id1]
  cdf[id2] <- log(1 - 0.5 * exp(-q[id2]))

  if (!log.p) cdf <- exp(cdf)
  cdf
}

qLA <- function(p) {
  qtf <- rep(NaN, length.out = length(p))

  # Positive density index
  id1 <- which(p <= 0.5 & p != 0)
  id2 <- which(p > 0.5 & p != 1)
  id3 <- which(p == 0)
  id4 <- which(p == 1)

  qtf[id1] <- log(2 * p[id1])
  qtf[id2] <- -log(2 - 2 * p[id2])
  qtf[id3] <- -Inf
  qtf[id4] <- Inf

  qtf
}
