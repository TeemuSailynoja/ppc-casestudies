library("reliabilitydiag")
library("ggplot2")
library("ggdist")
library("khroma")

plot_dotted_reliabilitydiag <- function(y,
                                        x,
                                        quantiles=100,
                                        region.method = "resampling",
                                        region.position = "diagonal",
                                        n.boot = 100,
                                        dot_scale = .25) {

  rdiag <- reliabilitydiag(y = y,
                           x = x,
                           region.method = region.method,
                           region.position = region.position,
                           n.boot = n.boot)

  ggplot() + autolayer(rdiag,
                       type = "miscalibration",
                       params_histogram = NA,
                       params_ribbon = list(
                         fill = colour("vibrant")(5)["cyan"],
                         alpha = .8
                       ),
                       params_CEPsegment = NA,
                       params_CEPline = list(
                         size = .5,
                         colour = colour("vibrant")(5)["orange"]
                       ),
                       params_CEPpoint = list(
                         size = 1,
                         colour = colour("vibrant")(5)["orange"]
                       )) +
    stat_dots(aes(x = x),
              color = "#666666",
              fill = "darkgrey",
               scale = dot_scale,
               quantiles = quantiles,
               inherit.aes = FALSE,
               alpha = .8) + ylab("CEP") + xlab("Forecast value")
}

plot_dotted_reliabilitydiag_old <- function(y,
                                        x,
                                        quantiles=100,
                                        region.method = "resampling",
                                        region.position = "diagonal",
                                        n.boot = 100,
                                        dot_scale = .25) {

  rdiag <- reliabilitydiag(y = y,
                           x = x,
                           region.method = region.method,
                           region.position = region.position,
                           n.boot = n.boot)

  p_rdiag <- autoplot(rdiag,
                      type = "miscalibration",
                      params_histogram = NA,
                      params_ribbon = list(
                        fill = colour("vibrant")(5)["blue"],
                        alpha = .25),
                      params_CEPsegment = NA,
                      params_CEPpoint = list(
                        size = 1,
                        colour = colour("vibrant")(5)["red"]
                        ),

                      )

  p_rdiag$layers <- c(stat_dots(aes(x = x),
                                color = "black",
                                fill = "white",
                                scale = dot_scale,
                                quantiles = quantiles,
                                inherit.aes = FALSE,
                                alpha = .8),
                      p_rdiag$layers)

  p_rdiag
}

geom_mark_suspicions <- function(p, llag = 1, llag_down = llag, ...) {
  xy <- layer_data(p)[, c("y", "x")]
  N <- nrow(xy)

  growth <- function(l) {head(diff(rowSums(xy), l) * N, N - l)}
  size <- function(l) {head(1 - rowSums(xy), N - l) * N}
  pr <- function(l) {(l / N) / (1 - head(xy$x, N - l))}

  filter_up <- head(rowSums(xy) < xy$x, N - llag) & (pbinom(q = growth(llag) - 1,
                                                            size = size(llag),
                                                            prob = pr(llag),
                                                            lower.tail = F,
                                                            log.p = T) < log(1 - 0.95 ** (1/N)))
  filter_up <- c(filter_up | c(rep(FALSE, llag), tail(filter_up, N - 2 * llag)), rep(FALSE, llag))

  group_up <- c(cumsum(pmax(0,diff(filter_up))) * filter_up[-(N - llag)], 0)

  filter_down <- head(rowSums(xy) < xy$x, N - llag_down) & (pbinom(q = growth(llag_down),
                                                                   size = size(llag_down),
                                                                   prob = pr(llag_down),
                                                                   log.p = T) < log(1 - 0.95 ** (1/N)))

  filter_down <- c(filter_down | c(rep(FALSE, llag_down), tail(filter_down, N - 2 * llag_down)), rep(FALSE, llag_down))

  group_down <- c(cumsum(pmax(0,diff(filter_down))) * filter_down[-(N - llag_down)], 0)

  group_down[group_down != 0] <- group_down[group_down != 0] + max(group_up)

  filter_over <- xy$y > layer_data(p, 2L)$y

  group_over <- filter_over + (filter_over * (xy$x > .1)) + (filter_over * (xy$x > .9))

  group_over[group_over != 0] <- group_over[group_over != 0] + max(group_down)

  filter_under <- xy$y < layer_data(p, 3L)$y

  group_under <- filter_under + (filter_under * (xy$x > .1)) + (filter_under * (xy$x > .9))

  group_under[group_under != 0] <- max(group_over) + group_under[group_under != 0]

  desc = c("",
           "local inflation of data",
           "local lack of data",
           "left tail thicker than expected",
           "mean smaller than expected",
           "right tail thinner than expected",
           "left tail thinner than expected",
           "mean larger than expected",
           "right tail thicker than expected")

  ggforce::geom_mark_rect(aes(filter = filter_up | filter_down | filter_over | filter_under,
                              group = group_up + group_down + group_over + group_under,
                              description = desc[1 + filter_up + 2 * filter_down +
                                                   (3 * filter_over + (filter_over * (xy$x > .1)) + (filter_over * (xy$x > .9))) +
                                                   (4 * filter_under + (filter_under * (xy$x > .1)) + (filter_under * (xy$x > .9)))],
                              colour = ), ...)
}

rooto_dots <- function(rooto,
                       shape = 21,
                       size = 3,
                       stroke = 1) {
  rooto$layers[[2]]$aes_params$linetype <- "dashed"
  rooto$layers[[2]]$aes_params$fill <-
    bayesplot::color_scheme_get()$mid
  rooto$layers[[2]]$aes_params$colour <-
    bayesplot::color_scheme_get()$mid_highlight
  rooto$layers <- rooto$layers[-1]
  rooto$data[, -1] <- rooto$data[, -1] ** 2

  rooto + ggplot2::geom_point(
    ggplot2::aes(
      x = xpos[ty != 0],
      y = ty[ty != 0],
      colour = ifelse((tylower <= ty) &
                        (ty <= tyupper), "Observed", FALSE),
      fill = ggplot2::after_scale(scales::alpha(colour, .5))
    ),
    data = subset(rooto$data, ty != 0),
    shape = shape,
    size = size,
    stroke = stroke
  ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Expected" = bayesplot::color_scheme_get()$mid,
        "Observed" = bayesplot::color_scheme_get()$dark,
        "FALSE" = "red"
      )
    ) +
    ggplot2::scale_y_continuous(
      trans = "sqrt",
      minor_breaks = function(limit) {
        floor(limit[1]):ceiling(limit[2])
      },
      breaks = function(limit) {
        round(seq(sqrt(limit[1]), ceiling(sqrt(limit[2])), length.out = 6) ^ 2)
      },
      expand = ggplot2::expansion(mult = c(0, .05))
    ) +
    ggplot2::ylab("Count") + ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_line(colour = "grey92"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
}

reflect_density <- function(dens, bounds, from, to) {
  # No adjustment is needed if no finite bounds are supplied
  if (all(is.infinite(bounds))) {
    return(dens)
  }

  # Estimate linearly with zero tails (crucial to account for infinite bound)
  f_dens <- stats::approxfun(
    x = dens$x, y = dens$y, method = "linear", yleft = 0, yright = 0
  )

  # Create a uniform x-grid inside `bounds`
  left <- max(from, bounds[1])
  right <- min(to, bounds[2])
  out_x <- seq(from = left, to = right, length.out = length(dens$x))

  # Update density estimation by adding reflected tails from outside `bounds`
  left_reflection <- f_dens(bounds[1] + (bounds[1] - out_x))
  right_reflection <- f_dens(bounds[2] + (bounds[2] - out_x))
  out_y <- f_dens(out_x) + left_reflection + right_reflection

  list(x = out_x, y = out_y)
}

pit_from_densityplot <- function(p, i, x, ggdist_layer = F) {
  require(sfsmisc)
  ld <- layer_data(p, i = i)
  if (ggdist_layer == TRUE) {
    return(unlist(lapply(x, function(x_i) if (x_i > min(x)) integrate.xy(ld$x, ld$thickness, b = x_i) else 0)))
  }
  unlist(lapply(x, function(x_i) if (x_i > min(x)) integrate.xy(ld$x, ld$y, b = x_i) else 0))
}
