library("shiny")
library("bayesplot")
library("ggplot2")
library("ggdist")
library("khroma")

good_theme <- bayesplot::theme_default(base_family = "Sans") + theme(
    axis.text = element_text(colour = "#666666", size = 12),
    axis.ticks = element_line(colour = "#666666"),
    title = element_text(colour = "#666666", size = 16),
    plot.subtitle = element_text(colour = "#666666", size = 14),
    legend.text = element_text(colour = "#666666", size = 12),
    legend.title = element_text(colour = "#666666", size = 14),
    axis.line = element_line(colour = "#666666"))

theme_set(good_theme)
bayesplot_theme_set(good_theme)
color_scheme_set(scheme = c(unname(colour("vibrant")(7)[c(3,2,5,4,6,1)])))

scale_colour_discrete = scale_colour_vibrant
scale_fill_discrete = scale_fill_vibrant

SEED <- 236543
set.seed(SEED)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Assessing KDE fit to data with ECDF plots"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("N",
                        "Number of observations:",
                        min = 100,
                        max = 500,
                        value = 250),
            sliderInput("n_na",
                        "Number of added values:",
                        min = 0,
                        max = 99,
                        value = 0),
            sliderInput("discontinuity_location",
                        "Location of added values:",
                        min = -1,
                        max = 3,
                        value = 0),
            sliderInput("n_rep",
                        "Number of posterior predictions:",
                        min = 0,
                        max = 10,
                        value = 5),
            width = 2
    ),
    mainPanel(
           plotOutput("densPlot"),
           plotOutput("ecdfPlot"),
           plotOutput("ecdfdiffPlot"),
           plotOutput("dotsPlot"),
           plotOutput("histPlot")
    )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    dataInput <- reactive({
        n_rep <- input$n_rep
        N <- input$N
        n_na <- input$n_na
        x <- rnorm(N, 1, .5)
        if (n_na > 0) {
            x[1:n_na] <- input$discontinuity_location
        }

        mu_0 <- 0
        nu <- 1
        alpha <- 1
        beta <- 1

        alpha_post <- alpha + N / 2
        beta_post <- beta + 0.5 * (N - 1) * var(x) +
            (N * nu * (mean(x) - mu_0) ** 2) / (2 * (nu + N))

        sd_post <- sqrt(1 / rgamma(n_rep, alpha_post, beta_post))
        nu_post <- nu + N
        mu_post <- rnorm(n_rep, nu * mu_0 + sum(x) / (nu * N), sd = sd_post / sqrt(nu_post) )

        x_post <- c(replicate(N, rnorm(n_rep, mu_post, sd_post)))
        rep_id <- rep(c(1:n_rep), each = N)
        p_dens <- ppc_dens_overlay(x, matrix(x_post, nrow = input$n_rep))
        fig_data <- layer_data(p_dens)
        KDEs <- data.frame(list(x = fig_data$x[fig_data$group == 1]))
        step_size <- KDEs$x[2] - KDEs$x[1]
        if (n_rep > 0) {
            for (id in 1:n_rep) {
                new = list()
                new[[paste("x_",id, sep="")]] = step_size * cumsum(fig_data[fig_data$group == id, ]$y)
                KDEs = cbind(KDEs, data.frame(new))
            }
        }
        KDEs = cbind(KDEs, data.frame(list(y = step_size * cumsum(layer_data(p_dens, 2L)[layer_data(p_dens, 2L)$group == 1,]$y))))
        list(x = x,
             x_post = x_post,
             ycol = layer_data(p_dens, 2L)$colour[1],
             yrepcol = layer_data(p_dens, 1L)$colour[1],
             KDEs = KDEs,
             linewidth = layer_data(p_dens)$linewidth[1],
             rep_id = rep_id
             )
    })

    output$densPlot <- renderPlot({
        data <- dataInput()
        ppc_dens_overlay(data$x, matrix(data$x_post, nrow = input$n_rep))
    })

    output$dotsPlot <- renderPlot({
        data <- dataInput()
        ggplot() +
            stat_dots(aes(x = data$x),
                      quantiles = 100,
                      fill = "transparent",
                      colour = data$ycol) +
            stat_density(aes(x = data$x),
                         geom = "line",
                         colour = data$ycol) +
            xlab("")
    })

    output$histPlot <- renderPlot({
        data <- dataInput()
        ggplot() +
            geom_histogram(aes(x = data$x),
                           colour = data$ycol,
                           fill = data$ycol,
                           alpha = .2)
    })

    output$ecdfPlot <- renderPlot({
        data <- dataInput()
        p_ecdf <- ppc_pit_ecdf(
            pit = sapply(data$x, function(x_i) data$KDEs$y[which.max(data$KDEs$x >= x_i)]),
            interpolate_adj = T)
        if (input$n_rep > 0) {
            for (idx in 1:input$n_rep) {
                gg_data = data.frame(
                    PIT = seq(0,1,length.out = length(data$x)),
                    ECDF = ecdf(
                        sapply(
                            data$x_post[data$rep_id == idx],
                            function(x_i) data$KDEs[, (1 + idx)][which.max(data$KDEs$x >= x_i)]))(seq(0,1,length.out = length(data$x)))
                )
                p_ecdf <- p_ecdf + geom_step(data = gg_data, aes(x = PIT, y = ECDF), colour = data$yrepcol, alpha = .7, linewidth = data$linewidth)
            }
            }
        p_ecdf
    }, bg = "transparent")

    output$ecdfdiffPlot <- renderPlot({
        data <- dataInput()
        p_ecdfd <- ppc_pit_ecdf(
            pit = sapply(data$x, function(x_i) data$KDEs$y[which.max(data$KDEs$x >= x_i)]),
            interpolate_adj = T, plot_diff = T)

        if (input$n_rep >0) {
            for (idx in 1:input$n_rep) {
                temp <- ecdf(
                    sapply(
                        data$x_post[data$rep_id == idx],
                        function(x_i) data$KDEs[ ,1 + idx][which.max(data$KDEs$x >= x_i)]))(seq(0,1,length.out = length(data$x)))
                gg_data = data.frame(
                    PIT = seq(0,1,length.out = length(data$x)),
                    ECDF = temp)
                p_ecdfd <- p_ecdfd + geom_step(data = gg_data, aes(x = PIT, y = ECDF - PIT), colour = data$yrepcol, alpha = .7, linewidth = data$linewidth)
            }
        }
        p_ecdfd
    })
}

# Run the application
shinyApp(ui = ui, server = server)
