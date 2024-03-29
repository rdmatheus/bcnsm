% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05b_mle.R
\name{bcnsm_fit}
\alias{bcnsm_fit}
\title{BCNSM Fit for Multivariate Positive Data}
\usage{
bcnsm_fit(
  y,
  association = c("unstructured", "uniform", "nonassociative"),
  copula = c("gaussian", "t", "slash", "hyp"),
  delta = NULL,
  margins = "bcno",
  control = control_fit(...),
  ...
)
}
\arguments{
\item{y}{a matrix with the sample observations.}

\item{association}{one of \code{"unstructured"} (default), \code{"uniform"}, or
\code{"nonassociative"}, which specify the association matrix of the model.}

\item{copula}{character; informs which normal scale mixture distribution should be used to generate
the NSM copula. Currently, the copulas available are: Gaussian (\code{"gaussian"}),
Student's \emph{t} (\code{"t"}), slash (\code{"slash"}), and hyperbolic (\code{"hyp"}).}

\item{delta}{possible extra parameter associated with the mixing distribution of the copula.
For example, the degrees of freedom of the \emph{t} copula, or the heavy tailness parameter of the
slash or the hyperbolic copula. To estimate the value of \code{delta} via profile log-likelihood,
use the function \code{\link{choose_copula}}.}

\item{margins}{a character or a character vector; specifies the marginal BCS distributions.
If all BCS margins are the same, it is sufficient to enter only one character. A table with
the current available BCS distributions can be seen in \code{\link{bcs}}.}

\item{control}{a list of control arguments specified via \code{\link{control_fit}}.}

\item{...}{further arguments passed to \code{\link{control_fit}}.}
}
\value{
The \code{bcnsm_fit} returns an object of class \code{"bcnsm"},
    which consists of a list with the following components:
\describe{
  \item{mu, sigma, lambda, nu}{vectors with the estimated values for \code{mu}, \code{sigma},
      \code{lambda}, and \code{nu}, respectively.}
  \item{gamma}{the estimated parameters of the association matrix, if any.}
  \item{margins}{a character vector with the marginal BCS distributions of the fit.}
  \item{association}{structure of the association matrix. It can be one of \code{"non-associative"},
      \code{"unstructured"}, or \code{"uniform"}.}
  \item{copula, delta}{\code{"copula"} is a character which informs which normal scale mixture distribution
      was used to generate the NSM copula and \code{"delta"} is the possible extra parameter associated with
      the copula.}
  \item{logLik}{log-likelihood of the fitted model.}
  \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector.
      Specifically, the inverse of the observed information matrix, obtained via numeric Hessian matrix.}
  \item{y}{the response matrix.}
  \item{optim_params}{control optimization parameters used by \code{\link{control_fit}}.}
  \item{nobs,d}{the number of observations in the sample and the dimension of the response variable, respectively.}
  \item{call}{ the function call.}
 }
}
\description{
Fit a multivariate Box-Cox symmetric distribution generated by an NSM copula via maximum
likelihood.
}
\examples{
\dontrun{
## Random sampling

### Sample size and dimension
n <- 200
d <- 4

### Association matrix
Gamma <- matrix(0.8, d, d)
diag(Gamma) <- 1

### Marginal specifications

margins <- c("bcno", "lt", "bct", "lno")

mu <- c(19, 20, 15, 20)
sigma <- c(0.2, 0.3, 0.4, 0.3)
lambda <- c(-1, NA, 1.6, NA)
nu <- c(NA, 4, 8, NA)

### Copula
copula <- "t"
delta <- 4

### Observations
set.seed(123) # Seed for reproducibility
y <- rbcnsm(n, mu, sigma, lambda, nu, Gamma, copula, delta, margins)
colnames(y) <- c("y1", "y2", "y3", "y4")

### Visualization (based on graphics::pairs functions)
mvplot(y)

## Fit with Gaussian copula and uniform structure
fit <- bcnsm_fit(y, association = "uniform", copula = "gaussian",
                 margins = c("bcno", "lt", "bct", "lno"))

class(fit)
methods(class = "bcnsm")

# Fit summaries
fit
summary(fit)

# Fit visualization

## Bivariate fit
plot(fit)

## Marginal fit
plot(fit, type = "margins")

## Transformed vectors
plot(fit, "epsilon")

# Choose the value of the extra parameter of the t copula (it can be slow)
fit_t <- choose_copula(fit, grid = 1:8, copula = "t")

## Final fit
final_fit <- fit_t[[4]]

final_fit
plot(final_fit)
plot(final_fit, type = "margins")
}
}
\references{
Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical properties and
    parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.

Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and applications to
    nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.

Medeiros, R. M. R. de, and Ferrari, S. L. P. (2023). Multivariate Box-Cox symmetric distributions
    generated by a normal scale mixture copula.
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
}
