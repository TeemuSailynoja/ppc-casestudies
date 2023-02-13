library("reliabilitydiag")
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
