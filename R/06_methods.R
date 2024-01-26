#' @name bcnsm-methods
#' @title Methods for \code{"bcnsm"} objects
#' @param x,object an object of class \code{bcnsm}.
#' @param digits number of digits in print methods.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

#  Variance-covariance matrix
#' @rdname bcnsm-methods
#' @param par character; specifies which submatrix of the asymptotic covariance matrix of the
#'     maximum likelihood estimators should be returned. The options are \code{"all"} (default),
#'     \code{"mu"}, \code{"sigma"}, \code{"lambda"}, \code{"nu"}, and \code{"gamma"}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
vcov.bcnsm <- function(object, par = c("all", "mu", "sigma", "lambda", "nu", "gamma"), ...) {

  par <- match.arg(par, c("all", "mu", "sigma", "lambda", "nu", "gamma"))
  covm <- object$vcov

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  d <- object$d

  names_mu <- paste0("mu", 1:d)
  names_sigma <- paste0("sigma", 1:d)
  names_lambda <- if (any(lambda_id)) paste0("lambda", (1:d)[lambda_id]) else NULL
  names_nu <- if (any(nu_id)) paste0("nu", (1:d)[nu_id]) else NULL

  if (object$association == "uniform") {

    names_gamma <- "gamma"

  } else {

    names_gamma <- vector()
    id <- which(upper.tri(diag(d)), arr.ind = TRUE)[order(which(upper.tri(diag(d)), arr.ind = TRUE)[, 1]), ]
    for (i in 1:nrow(id)) {
      names_gamma[i] <- paste0("gamma", id[i, 1], id[i, 2])
    }

  }

  colnames(covm) <- rownames(covm) <- c(names_mu, names_sigma, names_lambda, names_nu, names_gamma)

  par_id <- object$optim_params$par_id

  switch(par,
         "all" = covm,
         "mu" = covm[par_id$mu, par_id$mu],
         "sigma" = covm[par_id$sigma, par_id$sigma],
         "lambda" = covm[par_id$lambda, par_id$lambda],
         "nu" = covm[par_id$nu, par_id$nu],
         "gamma" = covm[par_id$gamma, par_id$gamma]
  )
}

# Log-likelihood
#' @rdname bcnsm-methods
#' @export
logLik.bcnsm <- function(object, ...) {

  structure(object$logLik,
            df = length(object$optim_params$par) + as.numeric(!is.null(object$delta)),
            class = "logLik")

}

# AIC
#' @export
#' @rdname bcnsm-methods
AIC.bcnsm <- function(object, ..., k = 2) {

  npars <- length(object$optim_params$par) + as.numeric(!is.null(object$delta))
  AIC <- -2 * object$logLik + k * npars

  class(AIC) <- "AIC"
  AIC

}

# Residuals
#' @name residuals.bcnsm
#' @title Extract Model Residuals for a BerG Regression
#'
#' @param object an \code{'bcnsm'} object.
#' @param type character; specifies which residual should be extracted. The available arguments are
#'     \code{"mahalanobis"} (default), and \code{"epsilon"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
residuals.bcnsm <- function(object, type = c("mahalanobis", "epsilon"), ...) {

  type <- match.arg(type, c("mahalanobis", "epsilon"))

  y <- object$y

  n <- object$nobs
  d <- object$d
  epsilon <- matrix(NA, n, d)
  margins <- object$margins
  copula <- object$copula
  delta <- object$delta

  copula_inf <- make_copula(copula, delta)
  qPSI <- copula_inf$qPSI

  mu <- object$mu
  sigma <- object$sigma
  lambda <- object$lambda
  nu <- object$nu

  for (j in 1:d) {

    epsilon[, j] <- qPSI(get(paste0("p", margins[j]))(q = y[, j],
                                                      mu = mu[j],
                                                      sigma = sigma[j],
                                                      lambda = lambda[j],
                                                      nu = nu[j]))

  }

  epsilon
  colnames(epsilon) <- colnames(y)

  Gamma <- get(object$association)(d)$Gamma(object$gamma)

  # Squared Mahalanobis distance
  mahalanobis <- Rfast::mahala(epsilon, rep(0L, d), Gamma)

  # Out
  res <- switch(type,
                "mahalanobis" = as.numeric(mahalanobis),
                "epsilon" = epsilon
  )

  res
}

# Print
#' @rdname bcnsm-methods
#' @export
print.bcnsm <- function(x, digits = 2L, ...) {

  y <- x$y

  n <- x$nobs
  d <- x$d
  margins <- x$margins
  association <- x$association
  gamma <- x$gamma
  copula <- x$copula
  delta <- x$delta

  # Parameters
  mu <- x$mu
  sigma <- x$sigma
  lambda <- x$lambda
  nu <- x$nu

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  # Names
  y_names <- if (is.null(colnames(y))) paste0("y", 1:d) else colnames(y)
  Margins <- toupper(margins)
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")

  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(crayon::cyan(
      "\nMultivariate Box-Cox Distribution with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    ))

    cat(crayon::cyan("Call:\n\n"))
    print(x$call)
    if (x$optim_params$convergence > 0) {

      cat("\nmodel did not converge\n")

    } else {

      cat(crayon::cyan("\nMarginal fit:\n\n"))

      comma_lambda <- comma_nu <- rep(",", d)
      lambda <- round(lambda, digits)
      nu <- round(nu, digits)

      lambda[!lambda_id] <- comma_lambda[!lambda_id] <- ""
      nu[!nu_id] <- comma_nu[!nu_id] <- ""

      for (j in 1:d) {

        cat(y_names[j], " ~ ", Margins[j], "(", round(mu[j], digits), ",",
            round(sigma[j], digits), comma_lambda[j], lambda[j], comma_nu[j], nu[j], ")\n",
            sep = ""
        )

      }

      if (length(gamma) > 0) {

        cat(crayon::cyan("\n", sub("(.)", "\\U\\1", x$association, perl = TRUE),
                         " association matrix:\n\n", sep = ""))

        Gamma <- get(x$association)(d)$Gamma(round(x$gamma, digits))
        Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
        colnames(Gamma) <- rownames(Gamma) <- y_names
        print(Gamma, na.print = ".")

      }

      cat(
        "\n---",
        crayon::cyan("\nlogLik:"), x$logLik, "|",
        crayon::cyan("AIC:"), stats::AIC(x), "|",
        crayon::cyan("BIC:"), stats::AIC(x, k = log(n)), "\n"
      )
    }

  } else {

    cat(
      "\nMultivariate Box-Cox Distribution with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    )

    cat("Call:\n\n")
    print(x$call)
    if (x$optim_params$convergence > 0) {

      cat("\nmodel did not converge\n")

    } else {

      cat("\nMarginal fit:\n\n")

      comma_lambda <- comma_nu <- rep(",", d)
      lambda <- round(lambda, digits)
      nu <- round(nu, digits)

      lambda[!lambda_id] <- comma_lambda[!lambda_id] <- ""
      nu[!nu_id] <- comma_nu[!nu_id] <- ""

      for (j in 1:d) {

        cat(y_names[j], ": ", Margins[j], "(", round(mu[j], digits), ",",
            round(sigma[j], digits), comma_lambda[j], lambda[j], comma_nu[j], nu[j], ")\n",
            sep = ""
        )

      }


      if (length(gamma) > 0) {

        cat("\n", sub("(.)", "\\U\\1", x$association, perl = TRUE), " association matrix:\n\n", sep = "")

        Gamma <- get(x$association)(d)$Gamma(round(x$gamma, digits))
        Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
        colnames(Gamma) <- rownames(Gamma) <- y_names
        print(Gamma, na.print = ".")
      }


      cat(
        "\n--------------------------------------------------",
        "\nlogLik:", x$logLik, "|",
        "AIC:", stats::AIC(x), "|",
        "BIC:", stats::AIC(x, k = log(n)), "\n"
      )
    }
  }



  invisible(x)
}

# Summary
#' @rdname bcnsm-methods
#' @export
summary.bcnsm <- function(object, digits = 4L, ...) {

  y <- object$y

  n <- object$nobs
  d <- object$d

  margins <- object$margins
  association <- object$association
  gamma <- object$gamma
  copula <- object$copula
  delta <- object$delta

  # Parameters
  mu <- object$mu
  sigma <- object$sigma
  lambda <- object$lambda
  nu <- object$nu

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  par_id <- object$optim_params$par_id

  # Names
  y_names <- if (is.null(colnames(y))) paste0("y", 1:d) else colnames(y)
  Margins <- toupper(margins)
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")

  # Summary for mu
  se_mu <- sqrt(diag(object$vcov)[par_id$mu])
  TAB_mu <- round(cbind(Estimate = c(mu), `Std. Error` = se_mu), digits)
  rownames(TAB_mu) <- colnames(y)

  # Summary for sigma
  se_sigma <- sqrt(diag(object$vcov)[par_id$sigma])
  TAB_sigma <- round(cbind(Estimate = c(sigma), `Std. Error` = se_sigma), digits)
  rownames(TAB_sigma) <- colnames(y)

  # Summary for lambda
  TAB_lambda <- NULL
  if (any(lambda_id)) {
    se_lambda <- sqrt(diag(object$vcov)[par_id$lambda])
    TAB_lambda <- round(cbind(Estimate = c(lambda[lambda_id]), `Std. Error` = se_lambda), digits)
    rownames(TAB_lambda) <- colnames(y)[lambda_id]
  }

  # Summary for nu
  TAB_nu <- NULL
  if (any(nu_id)) {
    se_nu <- sqrt(diag(object$vcov)[object$optim_params$par_id$nu])
    TAB_nu <- round(cbind(Estimate = object$nu[nu_id], `Std. Error` = se_nu), digits)
    rownames(TAB_nu) <- colnames(y)[nu_id]
  }

  gamma_id <- length(gamma) > 0

  # Summary for gamma
  TAB_gamma <- NULL
  if (gamma_id) {

    se_gamma <- sqrt(diag(object$vcov)[object$optim_params$par_id$gamma])
    TAB_gamma <- round(cbind(
      Estimate = gamma,
      `Std. Error` = se_gamma
      #`z value` = gamma / se_gamma,
      #`Pr(>|z|)` = 2 * stats::pnorm(abs(gamma / se_gamma), lower.tail = FALSE)
    ), digits)


    if (object$association == "uniform") {

      rownames(TAB_gamma) <- "gamma"

    } else {

      names_gamma <- vector()
      id <- which(upper.tri(diag(d)), arr.ind = TRUE)[order(which(upper.tri(diag(d)), arr.ind = TRUE)[, 1]), ]
      for (i in 1:nrow(id)) {
        names_gamma[i] <- paste0("gamma", id[i, 1], id[i, 2])
      }

      rownames(TAB_gamma) <- names_gamma

    }

  }

  npar <- length(object$optim_params$par)
  out <- list(
    mu = TAB_mu,
    sigma = TAB_sigma,
    lambda = TAB_lambda,
    nu = TAB_nu,
    gamma = TAB_gamma,
    margins = margins,
    association = association,
    copula = copula,
    y = y,
    d = d,
    delta = delta,
    logLik = object$logLik,
    AIC = stats::AIC(object),
    BIC = stats::AIC(object, k = log(n)),
    call = object$call
  )

  class(out) <- "summary.bcnsm"
  out
}

# Print summary
#' @rdname bcnsm-methods
#' @export
print.summary.bcnsm <- function(x, ...) {

  d <- x$d
  delta <- x$delta
  margins <- x$margins

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  # Names
  y_names <- if (is.null(colnames(x$y))) paste0("y", 1:d) else colnames(x$y)
  Margins <- toupper(x$margins)
  copula_name <- paste(toupper(substr(x$copula, 1, 1)), substr(x$copula, 2, nchar(x$copula)), sep = "")

  # Print
  if (requireNamespace("crayon", quietly = TRUE)) {

    cat(crayon::cyan(
      "\nMultivariate Box-Cox Distribution with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    ))

    cat(crayon::cyan("Call:\n\n"))
    print(x$call)

    cat(crayon::cyan("\nSummary for the marginal parameters:\n\n"))

    nu <- lambda <- matrix(NA, d, 2)
    colnames(nu) <- colnames(lambda) <- colnames(x$lambda)
    lambda[lambda_id, ] <- x$lambda
    nu[nu_id, ] <- x$nu

    for (j in 1:d) {

      cat(y_names[j], " ~ ", get(x$margins[j])()$name, " Distribution:\n", sep = "")

      TAB <- rbind(
        x$mu[j, ],
        x$sigma[j, ],
        if (lambda_id[j]) lambda[j, ],
        if (nu_id[j]) nu[j, ]
      )

      rownames(TAB) <- c("mu", "sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
      print(TAB)
      cat("\n\n")
    }

    if (!is.null(x$gamma)) {
      cat(crayon::cyan("Summary for the association parameters:\n\n", sep = ""))
      stats::printCoefmat(x$gamma)
    }

    article <- if (x$association == "non-associative") "a" else "an"
    cat(
      crayon::cyan("\nMarginal distributions:"), x$margins,
      crayon::cyan("\nDependence modeling:"), if (is.null(x$delta)) copula_name else paste0(copula_name, "(", x$delta, ")"),
      "copula with", article, tolower(x$association), "association matrix",
      crayon::cyan("\nLoglik:"), x$logLik, "|",
      crayon::cyan("AIC:"), x$AIC, "|",
      crayon::cyan("BIC:"), x$BIC
    )

  } else {

    cat(
      "\nMultivariate Box-Cox Distribution with",
      if (is.null(delta)) copula_name else paste0(copula_name, "(", round(delta, 2), ")"),
      "Copula\n\n"
    )

    cat("Call:\n\n")
    print(x$call)

    cat("\nSummary for the marginal parameters:\n\n")

    nu <- lambda <- matrix(NA, d, 4)
    colnames(nu) <- colnames(lambda) <- colnames(x$lambda)
    lambda[lambda_id, ] <- x$lambda
    nu[nu_id, 1:2] <- x$nu[, 1:2]

    for (j in 1:d) {

      cat(y_names[j], " ~ ", get(x$margins[j])()$name, " Distribution:\n", sep = "")

      TAB <- rbind(
        x$mu[j, ],
        x$sigma[j, ],
        if (lambda_id[j]) lambda[j, ],
        if (nu_id[j]) nu[j, ]
      )

      rownames(TAB) <- c("mu", "sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
      stats::printCoefmat(TAB, signif.legend = (j == d), na.print = "---")
      cat("\n\n")
    }

    if (!is.null(x$gamma)) {
      cat("Summary for the association parameters:\n\n", sep = "")
      stats::printCoefmat(x$gamma)
    }

    article <- if (x$association == "Non-associative") "a" else "an"
    cat(
      "\n---------------------------------------------------------------",
      "\nMarginal distributions:", x$margins,
      "\nDependence modeling:", if (is.null(x$delta)) copula_name else paste0(copula_name, "(", x$delta, ")"),
      "copula with", article, tolower(x$association), "association matrix",
      "\nLoglik:", x$logLik, "|",
      "AIC:", x$AIC, "|",
      "BIC:", x$BIC
    )
  }

  invisible(x)
}



globalVariables(c("theo", "emp", "marg", "grid_x", "grid_y", "prob", "density"))
# Plot
#' Visualization of the fit of the BCNSM distributions
#'
#'
#'
#' @param x an object of class \code{bcnsm}.
#' @param type character; specifies which graphical should be produced. The available options
#'     are \code{"response"} (default), \code{"margins"}, and \code{"epsilon"}.
#' @param outliers logical; used only when \code{type = "response"}. If \code{TRUE},
#'     possible outliers are highlighted in red.
#' @param alpha criterion according to the squared Mahalanobis distances that identifies a point as
#'     a possible outlier.
#' @param levels levels for contours plots.
#' @param panel A vector of the form \code{c(nr, nc)} with the number of rows and columns, where
#'     the figures will be drawn in an \code{nr-}by\code{-nc} array on the device.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
plot.bcnsm <- function(x, type = c("response", "margins", "epsilon"),
                       outliers = FALSE, alpha = 0.01,
                       panel = NULL,
                       levels = c(1, 1e-1, 1e-2, 1e-3, 1e-4), ...) {

  type <- match.arg(type, c("response", "margins", "epsilon"))

  y <- x$y
  epsilon <- stats::residuals(x, "epsilon")

  n <- x$nobs
  d <- x$d
  margins <- x$margins
  copula <- x$copula
  delta <- x$delta

  mu <- x$mu
  sigma <- x$sigma
  lambda <- x$lambda
  nu <- x$nu

  gamma <- x$gamma
  digits <- 4L
  if (is.null(panel)) panel <- c(ceiling(d / 4), 4)

  Gamma <- get(x$association)(d)$Gamma(round(x$gamma, digits))

  md <- Rfast::mahala(epsilon, rep(0L, d), Gamma)
  id_md <- get(paste0("maha_", copula))(md, delta, d) > 1 - alpha

  y_names <- if (is.null(colnames(y))) paste0("y", 1:d) else colnames(y)

  copula_inf <- make_copula(copula, delta)
  dmv <- get(paste0("dmv_", copula))
  dPSI <- copula_inf$dPSI
  pPSI <- copula_inf$pPSI
  qPSI <- copula_inf$qPSI

  op <- graphics::par()

  y_aux <- y
  y_names <- colnames(y)

  # With ggplot ------------------------------------------------------------------------------------
  if (requireNamespace("ggplot2", quietly = TRUE) &
      requireNamespace("GGally", quietly = TRUE)) {


    ## Response contour plot -----------------------------------------------------------------------
    if (type == "response") {

      ### Diagonal plots
      diag_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        xid <- which(colSums(y_aux - x) == 0)

        dBCS <- function(x) {
          get(paste0("d", margins[xid]))(x,
                                         mu = mu[xid],
                                         sigma = sigma[xid],
                                         lambda = lambda[xid],
                                         nu = nu[xid])
        }

        ggplot2::ggplot(data, mapping) +
          ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                  colour = 1, fill = "white",
                                  bins = ceiling(1 + 3.33 * log(n))) +
          ggplot2::geom_function(fun = dBCS, col = "#00AFBB")


      }

      ### Upper plots
      upper_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        i <- which(colSums(y_aux - x) == 0)
        j <- which(colSums(y_aux - y) == 0)

        corr <- Gamma[i, j]

        colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
        fill <- colFn(1000)[findInterval(corr, seq(-1, 1, length = 1000))]

        GGally::ggally_text(
          label = as.character(round(corr, 4)),
          mapping = ggplot2::aes(),
          xP = 0.5,
          yP = 0.5,
          cex = 2.5,
          color = 'black', ...) +
          ggplot2::theme_void() +
          ggplot2::theme(panel.background = ggplot2::element_rect(fill = fill))
      }

      ### Lower plots
      lower_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        i <- which(colSums(y_aux - x) == 0)
        j <- which(colSums(y_aux - y) == 0)

        Gamma_aux <- matrix(c(1, Gamma[i, j], Gamma[i, j], 1), 2, 2)

        mu_aux <- mu[c(i, j)]
        sigma_aux <- sigma[c(i, j)]
        lambda_aux <- lambda[c(i, j)]
        nu_aux <- nu[c(i, j)]

        grid <- expand.grid(seq(min(x) - 10, max(x) + 10, length.out = 200),
                            seq(min(y) - 10, max(y) + 10, length.out = 200))

        data_aux <- data.frame(grid_x = grid[, 1],
                               grid_y = grid[, 2],
                               prob = dbcnsm(cbind(grid[, 1], grid[, 2]),
                                             mu_aux, sigma_aux,
                                             lambda_aux, nu_aux, Gamma_aux,
                                             delta = delta, copula = copula,
                                             margins = margins[c(i, j)]))

        data_aux2 <- data.frame(x = x[id_md], y = y[id_md])


        out <- ggplot2::ggplot() +
          ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_contour(data = data_aux, mapping =  ggplot2::aes(x = grid_x, y = grid_y, z = prob),
                                col = "#00AFBB", breaks = levels) +
          ggplot2::labs(x = "", y = "")

        if (outliers) {
          out + ggplot2::geom_point(data = data_aux2, mapping = ggplot2::aes(x = x, y = y), col = "red")
        } else {
          out
        }



      }

      colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
      lab <- ggplot2::ggplot(data.frame(x = stats::runif(200, -1, 1), y = stats::runif(200, -1, 1),
                                        z = stats::runif(200, -1, 1)), ggplot2::aes(x, y, colour = z)) +
        ggplot2::geom_point() +
        ggplot2::scale_colour_gradient2("Fitted association \nparameter",
                                        low = colFn(200)[1], high = colFn(200)[200]) +
        ggplot2::theme(legend.title.align = 0.5,
                       legend.position = "top",
                       legend.key.height = ggplot2::unit(0.3, 'cm'),
                       legend.key.width = ggplot2::unit(1.5, 'cm'))

      GGally::ggpairs(as.data.frame(y), #ggplot2::aes(colour = gender),
                      upper = list(continuous = GGally::wrap(upper_func)),
                      lower = list(continuous = GGally::wrap(lower_func)),
                      diag = list(continuous = GGally::wrap(diag_func)),
                      legend = GGally::grab_legend(lab),
                      progress = FALSE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::theme(legend.position = "top")


    } else if (type == "margins") {

      ## Marginal distributions ------------------------------------------------------------------------
      theo_CDF <- emp_CDF <- matrix(NA, n, ncol(y))
      ks <- vector()
      for(j in 1:d){

        theo_CDF[ , j] <- get(paste0("p", margins[j]))(q = y[, j],
                                                       mu = mu[j],
                                                       sigma = sigma[j],
                                                       lambda = lambda[j],
                                                       nu = nu[j])

        ks[j] <- round(stats::ks.test(theo_CDF[, j], "punif")$p.value, 2)

        emp_CDF[, j] <- stats::ecdf(y[, j])(y[, j])

        id <- order(theo_CDF[, j])

        theo_CDF[, j] <- theo_CDF[id, j]
        emp_CDF[, j] <- emp_CDF[id, j]

      }

      dKS <- sfsmisc::KSd(n)
      aux_data <- data.frame(theo = c(theo_CDF),
                             emp = c(emp_CDF),
                             lower = c(emp_CDF - dKS),
                             upper = c(emp_CDF + dKS),
                             marg = factor(rep(colnames(y), each = n),
                                           levels = colnames(y)),
                             ks = rep(paste0("KS p-value: ", ks), each = n))

      positions <- data.frame(x = c(rbind(-1, theo_CDF, 2, apply(theo_CDF, 2, rev))),
                              y = c(rbind(-1 - dKS, theo_CDF - dKS, 2 + dKS, apply(theo_CDF, 2, rev) + dKS)),
                              marg = factor(rep(colnames(y), each = 2 * n + 2),
                                            levels = colnames(y)))

      ggplot2::ggplot() +
        ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = marg), positions, fill = "#cceff1") +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::geom_point(ggplot2::aes(x = theo, y = emp, group = marg), aux_data) +
        ggplot2::geom_abline(intercept = 0, slope = 1, col = "#00AFBB", lty = 1) +
        ggplot2::geom_text(ggplot2::aes(x = 0.8, y = 0, label = ks, group = marg), aux_data, size = 3) +
        ggplot2::facet_wrap(~ marg, nrow = panel[1], ncol = panel[2]) +
        ggplot2::labs(x = "Fitted distribution function",
                      y = "Empirical distribution function", size = 1)


    } else if (type == "epsilon") {
      ## Epsilon contour plot ----------------------------------------------------------------------


      ### Diagonal plots
      diag_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)

        ggplot2::ggplot(data, mapping) +
          ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                  colour = 1, fill = "white",
                                  bins = ceiling(1 + 3.33 * log(n))) +
          ggplot2::geom_function(fun = dPSI, col = "#00AFBB")


      }

      ### Upper plots
      upper_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        i <- which(colSums(epsilon - x) == 0)
        j <- which(colSums(epsilon - y) == 0)

        corr <- Gamma[i, j]

        colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
        fill <- colFn(1000)[findInterval(corr, seq(-1, 1, length = 1000))]

        GGally::ggally_text(
          label = as.character(round(corr, 4)),
          mapping = ggplot2::aes(),
          xP = 0.5,
          yP = 0.5,
          cex = 2.5,
          color = 'black', ...) +
          ggplot2::theme_void() +
          ggplot2::theme(panel.background = ggplot2::element_rect(fill = fill))
      }

      ### Lower plots
      lower_func <- function(data, mapping, ...){

        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        i <- which(colSums(epsilon - x) == 0)
        j <- which(colSums(epsilon - y) == 0)

        Gamma_aux <- matrix(c(1, Gamma[i, j], Gamma[i, j], 1), 2, 2)

        mu_aux <- mu[c(i, j)]
        sigma_aux <- sigma[c(i, j)]
        lambda_aux <- lambda[c(i, j)]
        nu_aux <- nu[c(i, j)]

        grid <- expand.grid(seq(min(x) - 10, max(x) + 10, length.out = 200),
                            seq(min(y) - 10, max(y) + 10, length.out = 200))

        data_aux <- data.frame(grid_x = grid[, 1],
                               grid_y = grid[, 2],
                               prob = dmv(cbind(grid[, 1], grid[, 2]),
                                          mu = rep(0, 2L),
                                          Sigma = Gamma_aux,
                                          delta = delta))


        ggplot2::ggplot() +
          ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_contour(data = data_aux, mapping =  ggplot2::aes(x = grid_x, y = grid_y, z = prob),
                                col = "#00AFBB", breaks = levels) +
          ggplot2::labs(x = "", y = "")


      }

      colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
      lab <- ggplot2::ggplot(data.frame(x = stats::runif(200, -1, 1), y = stats::runif(200, -1, 1),
                                        z = stats::runif(200, -1, 1)), ggplot2::aes(x, y, colour = z)) +
        ggplot2::geom_point() +
        ggplot2::scale_colour_gradient2("Fitted association \nparameter",
                                        low = colFn(200)[1], high = colFn(200)[200]) +
        ggplot2::theme(legend.title.align = 0.5,
                       legend.position = "top",
                       legend.key.height = ggplot2::unit(0.3, 'cm'),
                       legend.key.width = ggplot2::unit(1.5, 'cm'))

      GGally::ggpairs(as.data.frame(epsilon), #ggplot2::aes(colour = gender),
                      upper = list(continuous = GGally::wrap(upper_func)),
                      lower = list(continuous = GGally::wrap(lower_func)),
                      diag = list(continuous = GGally::wrap(diag_func)),
                      legend = GGally::grab_legend(lab),
                      progress = FALSE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::theme(legend.position = "top")


    }

    # With R base plots ------------------------------------------------------------------------------
  } else {

    if (type == "response") {

      ## Response contour plot -----------------------------------------------------------------------

      ## put histograms on the diagonal
      panel.hist <- function(x, ...) {

        xid <- which(colSums(y_aux - x) == 0)


        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))
        h <- graphics::hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$density
        graphics::par(usr = c(usr[1:2], 0, max(y) + 0.5 * diff(range(y))))
        graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)

        dBCS <- function(x) {
          get(paste0("d", margins[xid]))(x, mu = mu[xid], sigma = sigma[xid],
                                         lambda = lambda[xid], nu = nu[xid])
        }

        graphics::curve(dBCS(x), add = TRUE, col = "#00AFBB", lwd = 2, ...)

      }


      ## put an association measure on the upper panels
      panel.cor <- function(x, y, ...) {

        xid <- which(colSums(y_aux - x) == 0)
        yid <- which(colSums(y_aux - y) == 0)

        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))

        r <- Gamma[xid, yid]

        colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate = "spline")
        fill <- colFn(1000)[findInterval(r, seq(-1, 1, length = 1000))]

        graphics::par(usr = c(0, 1, 0, 1))
        txt <- format(c(r, 0.123456789), digits = 2)[1]
        txt <- paste0(txt)
        graphics::rect(0, 0, 1, 1, col = fill, border = "black")
        graphics::text(0.5, 0.5, txt, cex = 1.2)

      }

      ## put scatterplots with countour lines on the lower panels
      panel.scatter <- function(x, y, ...) {

        xid <- which(colSums(y_aux - x) == 0)
        yid <- which(colSums(y_aux - y) == 0)

        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))
        graphics::points(x, y, pch = 16)

        mu_aux <- mu[c(xid, yid)]
        sigma_aux <- sigma[c(xid, yid)]
        lambda_aux <- lambda[c(xid, yid)]
        nu_aux <- nu[c(xid, yid)]
        Gamma_aux <- matrix(c(1, Gamma[xid, yid], Gamma[xid, yid], 1), 2, 2)

        faux <- function(x, y) {

          dbcnsm(
            cbind(x, y), mu_aux, sigma_aux, lambda_aux, nu_aux, Gamma_aux,
            copula, delta, margins[c(xid, yid)]
          )

        }

        z <- outer(seq(min(x), max(x), length.out = 200),
                   seq(min(y), max(y), length.out = 200), faux)

        graphics::contour(seq(min(x), max(x), length.out = 200),
                          seq(min(y), max(y), length.out = 200), z,
                          levels = levels, col = "#00AFBB", lwd = 2, add = TRUE
        )



      }


      graphics::pairs(y, lower.panel = panel.scatter,
                      upper.panel = panel.cor,
                      diag.panel = panel.hist,
                      labels = y_names, cex.labels = 1.1, font.labels = 2)




    } else if (type == "margins") {

      ## Marginal distributions --------------------------------------------------------------------

      graphics::par(mfrow = panel)

      # P-P plot
      for (j in 1:d) {

        u <- stats::ecdf(y[, j])(y[, j])
        v <- get(paste0("p", margins[j]))(q = y[, j], mu = mu[j], sigma = sigma[j],
                                          lambda = lambda[j], nu = nu[j])

        dKS <- sfsmisc::KSd(n)
        id <- order(v)

        plot(v[id], u[id], main = y_names[j], pch = 16, xlim = c(0, 1), ylim = c(0, 1),
             xlab = "Fitted distribution function", ylab = "Empirical distribution function", type = "n")
        graphics::polygon(x = c(-1, u[id], 2, rev(u[id])),
                          y = c(-1 - dKS, u[id] - dKS, 2 + dKS, rev(u[id]) + dKS),
                          col = "#cceff1", border = NA)
        graphics::points(v, u, pch = 16)
        graphics::abline(0, 1, col = "#00AFBB", lwd = 2)
        graphics::box()
        graphics::legend("bottomright",
                         paste("KS p-value:", suppressWarnings(round(stats::ks.test(v, "punif", 0, 1)$p.value, 2))),
                         bty = "n", cex = 0.75)

      }

      graphics::par(mfrow = op$mfrow)

    } else if (type == "epsilon") {

      ## Epsilon contour plot ----------------------------------------------------------------------------

      ## put histograms on the diagonal
      panel.hist <- function(x, ...) {

        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))
        graphics::par(usr = c(usr[1:2], 0, 0.6))
        h <- graphics::hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$density
        graphics::par(usr = c(usr[1:2], 0, max(max(y), dPSI(0L)) + 0.5 * diff(range(y))))
        graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
        graphics::curve(dPSI(x), add = TRUE, col = "#00AFBB", lwd = 2)

      }

      ## put kendall's tau on the upper panels
      panel.cor <- function(x, y, ...) {

        xid <- which(colSums(epsilon - x) == 0)
        yid <- which(colSums(epsilon - y) == 0)

        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))

        r <- Gamma[xid, yid]

        colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate = "spline")
        fill <- colFn(1000)[findInterval(r, seq(-1, 1, length = 1000))]

        graphics::par(usr = c(0, 1, 0, 1))
        txt <- format(c(r, 0.123456789), digits = 2)[1]
        txt <- paste0("Fitted\n Kendall's tau:\n", txt)
        graphics::rect(0, 0, 1, 1, col = fill, border = "black")
        graphics::text(0.5, 0.5, txt, cex = 1.2)
      }

      ## put scatterplots with countour lines on the lower panels
      panel.scatter <- function(x, y, ...) {

        xid <- which(colSums(epsilon - x) == 0)
        yid <- which(colSums(epsilon - y) == 0)

        faux <- function(x, y, Gamma) {
          if (copula == "gaussian") {
            dmv_gaussian(cbind(x, y), rep(0, 2), Gamma)
          } else if (copula == "t") {
            dmv_t(cbind(x, y), rep(0, 2), Gamma, delta)
          } else if (copula == "slash") {
            dmv_slash(cbind(x, y), rep(0, 2), Gamma, delta)
          } else {
            dmv_hyp(cbind(x, y), rep(0, 2), Rfast::cholesky(Gamma), delta)
          }

        }

        Gamma_aux <- matrix(c(1, Gamma[xid, yid], Gamma[xid, yid], 1), 2, 2)
        z <- outer(seq(min(x), max(x), length.out = 200),
                   seq(min(y), max(y), length.out = 200), faux,
                   Gamma = Gamma_aux)

        usr <- graphics::par("usr")
        on.exit(graphics::par(usr = usr))
        graphics::points(x, y, pch = 16)
        graphics::contour(seq(min(x), max(x), length.out = 200),
                          seq(min(y), max(y), length.out = 200), z,
                          levels = levels, col = "#00AFBB", lwd = 2, add = TRUE
        )
      }

      graphics::pairs(epsilon,
                      lower.panel = panel.cor,
                      upper.panel = panel.scatter,
                      diag.panel = panel.hist,
                      labels = y_names, cex.labels = 1.1, font.labels = 2
      )

    }



  }

}
