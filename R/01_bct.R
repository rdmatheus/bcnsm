# The Box-Cox t distribution -----------------------------------------------------------------------

#' @name bct
#'
#' @title The Box-Cox t Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the Box-Cox t distribution with parameters \code{mu},
#'     \code{sigma}, \code{lambda}, and \code{nu}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lambda vector of real-valued skewness parameters. If \code{lambda = 0}, the Box-Cox
#'     t distribution reduces to the log-t distribution with parameters
#'     \code{mu}, \code{sigma}, and \code{nu} (see \code{\link{lt}}).
#' @param nu strictly positive heavy-tailedness parameter.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'
#' @return \code{dbct} returns the density, \code{pbct} gives the distribution function,
#'     \code{qbct} gives the quantile function, and \code{rbct} generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}, with an warning.
#'
#'     The length of the result is determined by \code{n} for \code{rbct}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Rigby, R. A., Stasinopoulos, D.M. (2006). Using the Box-Cox t
#'     distribution in GAMLSS to model skewness and kurtosis. \emph{Statistical Model}, 6, 209-229
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @references Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'     applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' mu <- 8
#' sigma <- 1
#' lambda <- 2
#' nu <- 4
#'
#' # Sample generation
#' x <- rbct(10000, mu, sigma, lambda, nu)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Box-Cox t Distribution", col = "white")
#' curve(dbct(x, mu, sigma, lambda, nu), add = TRUE, col = 2, lwd = 2)
#' legend("topleft", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Box-Cox t Distribution", ylab = "Distribution function")
#' curve(pbct(x, mu, sigma, lambda, nu), add = TRUE, col = 2, lwd = 2)
#' legend("topleft", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Box-Cox t Distribution"
#' )
#' curve(qbct(x, mu, sigma, lambda, nu), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL

# Density
#' @rdname bct
#' @export
dbct <- function(x, mu, sigma, lambda, nu, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id1 <- which(x > 0 & lambda != 0 & !is.nan(pmf))
  id2 <- which(x > 0 & lambda == 0 & !is.nan(pmf))

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  z[id1] <- ((x[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(x[id2] / mu[id2]) / sigma[id2]


  pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + stats::dt(z[id1], df = nu[id1], log = TRUE) -
    stats::pt(1 / (sigma[id1] * abs(lambda[id1])), df = nu[id1], log.p = TRUE) -
    lambda[id1] * log(mu[id1]) - log(sigma[id1])

  pmf[id2] <- stats::dt(z[id2], df = nu[id2], log = TRUE) - log(sigma[id2] * x[id2])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname bct
#' @export
pbct <- function(q, mu, sigma, lambda, nu, lower.tail = TRUE) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda != 0 & nu > 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0 & nu > 0)

  z[id1] <- ((q[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(q[id2] / mu[id2]) / sigma[id2]

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda <= 0 & nu > 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0 & nu > 0)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id1] <- stats::pt(z[id1], df = nu[id1]) / stats::pt(1 / (sigma[id1] * abs(lambda[id1])), df = nu[id1])
  cdf[id2] <- (stats::pt(z[id2], df = nu[id2]) - stats::pt(-1 / (sigma[id2] * lambda[id2]), df = nu[id2])) /
    stats::pt(1 / (sigma[id2] * lambda[id2]), df = nu[id2])

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname bct
#' @export
qbct <- function(p, mu, sigma, lambda, nu, lower.tail = TRUE) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(lambda), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda <= 0 & nu > 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda > 0 & nu > 0)

  zp[id1] <- stats::qt(p[id1] * stats::pt(1 / (sigma[id1] * abs(lambda[id1])), df = nu[id1]), df = nu[id1])
  zp[id2] <- stats::qt(1 - (1 - p[id2]) * stats::pt(1 / (sigma[id2] * abs(lambda[id2])), df = nu[id2]), df = nu[id2])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0 & nu > 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0 & nu > 0)
  id3 <- which(p == 0 & mu > 0 & sigma > 0 & nu > 0)
  id4 <- which(p == 1 & mu > 0 & sigma > 0 & nu > 0)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * zp[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * zp[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname bct
#' @export
rbct <- function(n, mu, sigma, lambda, nu) {
  u <- stats::runif(n)
  qbct(u, mu, sigma, lambda, nu)
}

# BCS class
bct <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "bct"

  # Name
  out$name <- "Box-Cox t"

  # Number of parameters
  out$npar <- 4

  # Extra parameter
  out$extrap <- TRUE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCT(mu.link = "log"), trace = FALSE), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")),
        stats::coef(gamlss_fit, "nu"), min(exp(stats::coef(gamlss_fit, "tau")), 20))
    } else {

      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)

      mu0 <- stats::median(x)
      sigma0 <- asinh(CV / 1.5) * stats::qlogis(0.75)

      z <- log(x / mu0) / sigma0

      grid <- seq(1, 20, 1)
      upsilon <- function(nu){
        cdf <- sort(stats::dt(z, nu))
        temp <- stats::qqnorm(stats::qnorm(cdf), plot.it = FALSE)
        mean(abs(sort(temp$x) - sort(temp$y)))
      }

      out <- apply(matrix(grid), 1, upsilon)
      nu0 <- grid[which.min(out)]

      c(mu0, sigma0, 0L, nu0)
    }

  }

  structure(out, class = "bcs")
}


# The Log-t distribution -------------------------------------------------------------------

#' @name lt
#'
#' @title The Log-t Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the log-t distribution with parameters \code{mu},
#'     \code{sigma}, and \code{nu}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param nu strictly positive heavy-tailedness parameter.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... further arguments.
#'
#' @details A random variable X has a log-t distribution with parameter \code{mu} and
#'     \code{sigma} if log(X) follows a Student's t distribution with location parameter \code{log(mu)}
#'     and dispersion parameter \code{sigma}. It can be showed that \code{mu} is the median of X.
#'
#' @return \code{dlt} returns the density, \code{plt} gives the distribution
#'     function, \code{qlt} gives the quantile function, and \code{rlt}
#'     generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}.
#'
#'     The length of the result is determined by \code{n} for \code{rlt}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' mu <- 8
#' sigma <- 1
#' nu <- 4
#'
#' # Sample generation
#' x <- rlt(10000, mu, sigma, nu)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Log-t Distribution", col = "white")
#' curve(dlt(x, mu, sigma, nu), add = TRUE, col = 2, lwd = 2)
#' legend("topright", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Log-t Distribution", ylab = "Distribution function")
#' curve(plt(x, mu, sigma, nu), add = TRUE, col = 2, lwd = 2)
#' legend("bottomright", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Log-t Distribution"
#' )
#' curve(qlt(x, mu, sigma, nu), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL

# Density
#' @rdname lt
#' @export
dlt <- function(x, mu, sigma, nu, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(nu)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0 | nu <= 0)] <- NaN

  # Positive density index
  id <- which(x > 0 & !is.nan(pmf))

  # Transformations
  z <- rep(NaN, length.out = maxl)

  z[id] <- log(x[id] / mu[id]) / sigma[id]


  pmf[id] <- stats::dt(z[id], df = nu[id], log = TRUE) - log(sigma[id] * x[id])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname lt
#' @export
plt <- function(q, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(nu)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id <- which(q > 0 & mu > 0 & sigma > 0 & nu > 0)

  z[id] <- log(q[id] / mu[id]) / sigma[id]

  cdf <- rep(NaN, length.out = maxl)
  cdf[id] <- stats::pt(z[id], df = nu[id])

  cdf[which(q <= 0 & mu > 0 & sigma > 0 & nu > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname lt
#' @export
qlt <- function(p, mu, sigma, nu, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(nu)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  nu <- rep(nu, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & nu > 0)

  zp[id] <- stats::qt(p[id], df = nu[id])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & nu > 0)
  id2 <- which(p == 0 & mu > 0 & sigma > 0 & nu > 0)
  id3 <- which(p == 1 & mu > 0 & sigma > 0 & nu > 0)

  qtf[id1] <- exp(log(mu[id1]) + sigma[id1] * zp[id1])
  qtf[id2] <- 0
  qtf[id3] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname lt
#' @export
rlt <- function(n, mu, sigma, nu) {
  exp(log(mu) + sigma * stats::rt(n, df = nu))
}

# BCS class
lt <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "lt"

  # Name
  out$name <- "Log-t"

  # Number of parameters
  out$npar <- 3

  # Extra parameter
  out$extrap <- TRUE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCT(mu.link = "log"),
                                                      trace = FALSE, nu.fix = TRUE, nu.start = 0L), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")),
        min(exp(stats::coef(gamlss_fit, "tau")), 20))
    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)

      mu0 <- stats::median(x)
      sigma0 <- asinh(CV / 1.5) * stats::qlogis(0.75)

      z <- log(x / mu0) / sigma0

      grid <- seq(1, 20, 1)
      upsilon <- function(nu){
        cdf <- sort(stats::dt(z, nu))
        temp <- stats::qqnorm(stats::qnorm(cdf), plot.it = FALSE)
        mean(abs(sort(temp$x) - sort(temp$y)))
      }

      out <- apply(matrix(grid), 1, upsilon)
      nu0 <- grid[which.min(out)]

      c(mu0, sigma0, nu0)

    }

  }



  structure(out, class = "bcs")
}
