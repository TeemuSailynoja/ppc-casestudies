library(ggplot2)
library(patchwork)
library(ggtext)

theme_set(theme_minimal(base_size = 14))

#' Adjust coverage parameter for single sample using the optimization method.
#' @param N Length of sample.
#' @param K Number of equally spaced evaluation points (1:K / K). Defaults to N.
#' @param prob Desired simultaneous coverage (0,1).
#' @return The adjusted coverage parameter yielding the desired simultaneous
#'  coverage of the ECDF traces.
#' @noRd
adjust_gamma_optimize <- function(N, K, prob) {
  target <- function(gamma, prob, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0, z)
    z2 <- c(z, 1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], 1)

    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probabilities known
    # to be 0 and 1 for the starting value z1 = 0.
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      p_int <- p_interior(
        p_int = p_int,
        x1 = x1,
        x2 = x2_lower[i]: x2_upper[i],
        z1 = z1[i],
        z2 = z2[i],
        N = N
      )
      x1 <- x2_lower[i]:x2_upper[i]
    }
    return(abs(prob - sum(p_int)))
  }
  optimize(target, c(0, 1 - prob), prob = prob, N = N, K = K)$minimum
}

#' A helper function for 'adjust_gamma_optimize' defining the probability that
#' a scaled ECDF stays within the supplied bounds between two evaluation points.
#' @param p_int For each value in x1, the probability that the ECDF has stayed
#' within the bounds until z1 and takes the value in x1 at z1.
#' @param x1 Vector of scaled ECDF values at the left end of the interval, z1.
#' @param x2 Vector of scaled ECDF values at the right end of the interval, z2.
#' @param z1 Left evaluation point in [0,1]
#' @param z2 Right evaluation point in [0,1] with z2 > z1.
#' @param N Total number of values in the sample.
#' @return A vector containing the probability to transitioning from the values
#' in x1 to each of the values in x2 weighted by the probabilities in p_int.
#' @noRd
p_interior <- function(p_int, x1, x2, z1, z2, N) {
  # Ratio between the length of the evaluation interval and the total length of
  # the interval left to cover by ECDF.
  z_tilde <- (z2 - z1) / (1 - z1)
  # Number of samples left to cover by ECDF.
  N_tilde <- rep(N - x1, each = length(x2))

  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  # Pobability of each transition from a value in x1 to a value in x2.
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)
  rowSums(p_x2_int)
}

#' Compute simultaneous confidence intervals with the given adjusted coverage
#'  parameter gamma.
#' @param gamma Adjusted coverage parameter for the marginal distribution
#'  (binomial for PIT values and hypergeometric for rank transformed chains).
#' @param N Sample length.
#' @param K Number of uniformly spaced evaluation points.
#' @return A list with upper and lower confidence interval values at the
#' evaluation points.
#' @noRd
ecdf_confidence_intervals <- function(gamma, N, K) {
  lims <- list()
  z <- seq(0, 1, length.out = K + 1)
  lims$lower <- qbinom(gamma / 2, N, z)
  lims$upper <- qbinom(1 - gamma / 2, N, z)
  lims
}

K = 200 # Number of evaluation points for ECDF
N = 500 # Sample size
M = 100 # number of simulations
prob = 0.95

# Limits for detecting non-uniformity
lims_outer <- ecdf_confidence_intervals(
  (gamma_outer <- adjust_gamma_optimize(N, K, prob + .5 * (1 - prob))),
  N = N, K = K)

# Limits for detecting too-uniformity
lims_inner <- ecdf_confidence_intervals(
  (gamma_inner <- adjust_gamma_optimize(N, K, .5 * (1 - prob))),
  N = N, K = K)

# Test if a given sample is not uniform,.
not_unif <- function(x, lims) {
  any(
    head(tail(x > lims$upper,-1),-1)| head(tail(x < lims$lower, -1),-1)
    )
}

too_unif <- function(x, lims) {
  all(
    head(tail(x <= lims$upper,-1),-1) & head(tail(x >= lims$lower, -1),-1)
    )
}

# Check that rejection rate is close to 'prob'
mean(replicate(M, {
  X <- colSums(outer(runif(N), (0:K / K), "<="))
  not_unif(X, lims_outer) | too_unif(X, lims_inner)
  }
  ))


#

p <- ggplot(
  data = data.frame(
    ECDF = c(c(sapply(lims_outer, c)),
          c(sapply(lims_inner, c))) / N,
    PIT = rep(0:K / K, 4),
    colour = "lims"
  )) +
    aes(x = PIT, y = ECDF - PIT, colour = colour) +
    geom_step(aes(group = rep(
      c("upper_outer", "lower_outer", "upper_inner", "lower_inner"),
      each = K + 1)), show.legend = F) +
  scale_colour_manual(
        name="",
        values = c("unif" = "grey",
                   "not_unif"= "red",
                   "too_unif"= "blue"),
        labels = c("unif" = "Uniform",
                   "not_unif"= "Not Uniform",
                   "too_unif"= "Too Uniform"),
        drop = FALSE) +
  theme(legend.position = "none",
        plot.title.position = "plot")
p_of <- p + labs(subtitle = "Half bandwitdth")
p_uf <- p + labs(subtitle = "Double bandwitdth")
p <- p + labs(title = "**KDE fit to data**<br>
              <span style='font-size:14pt;'>Below, we sample observations from the standard normal distribution, and then
              evaluate the ECDF of a KDE to that sample at the observed values.
              If the KDE represents the true distribution, the ECDF values should be <span style='color:#808080;'>**UNIFORMLY DISTRIBUTED**</span>.
              If the KDE is underfitting to the data, the ECDF is often <span style='color:#FF0000;'>**NOT UNIFORM**</span>.
              An overfit to the data causes te ECDF to be <span style='color:#0000FF;'>**TOO UNIFORM**</span>. ",
              subtitle = "Default bandwitdth") +
  theme(
    plot.title = element_textbox_simple(
      size = 17, lineheight = 1, hjust = 0, padding = margin(0, 0, 10, 0),
      width = grid::unit(.5, "npc"),
      ))

kernel <- "gaussian"

for (m in 1:M) {
  # simulate a sample
  x <- rnorm(N)

  # Fit thee KDEs (default, half bw, and double bw)
  f <- density(x,
               bw = "SJ", kernel = kernel)
  f_of <- density(x, kernel = kernel,
                  bw = "SJ", adjust = .5)
  f_uf <- density(x, kernel = kernel,
                  bw = "SJ", adjust = 2)

  # d_x for each of the KDEs. Used below in PIT calculation.
  delta <- f$x[2] - f$x[1]
  delta_of <- f_of$x[2] - f_of$x[1]
  delta_uf <- f_uf$x[2] - f_uf$x[1]

  # Compute PIT of observations
  F_x   <- ecdf(sapply(x, function(x_i) delta * sum(f$y[f$x < x_i])))(0:K / K)
  Fof_x <- ecdf(sapply(x, function(x_i) delta_of * sum(f_of$y[f_of$x < x_i])))(0:K / K)
  Fuf_x <- ecdf(sapply(x, function(x_i) delta_uf * sum(f_uf$y[f_uf$x < x_i])))(0:K / K)

  # Add each ECDF difference plot to the subplots.
  p <- p + geom_step(
    data = data.frame(
            ECDF = F_x,
            PIT = 0:K / K,
            colour = rep({if (not_unif(N * F_x, lims_outer)) "not_unif" else
              if (too_unif(N * F_x, lims_inner)) "too_unif" else
                "unif"}, K + 1)
            ),
    alpha = .3)

  p_of <- p_of + geom_step(
    data = data.frame(
      ECDF = Fof_x,
      PIT = 0:K / K,
      colour = rep({if (not_unif(N * Fof_x, lims_outer)) "not_unif" else
        if (too_unif(N * Fof_x, lims_inner)) "too_unif" else
          "unif"}, K + 1)
    ),
    alpha = .3)

  p_uf <- p_uf + geom_step(
    data = data.frame(
      ECDF = Fuf_x,
      PIT = 0:K / K,
      colour = rep({if (not_unif(N * Fuf_x, lims_outer)) "not_unif" else
        if (too_unif(N * Fuf_x, lims_inner)) "too_unif" else
          "unif"}, K + 1)
    ),
    alpha = .3)
}

# Show the resulting plot.
p / p_of / p_uf
