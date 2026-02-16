library(shiny)
library(shinyjs)
library(forecast)
library(ggplot2)
library(tseries)
library(urca)

ui <- fluidPage(
  useShinyjs(), 
  titlePanel("ARIMA Identification & Modelling"),
  
  # Full-screen placeholder text
  fluidRow(
    column(12, 
           wellPanel(
             style = "background-color: #f8f9fa; border-left: 5px solid #007bff;",
             tags$p("Iteratively identify ARIMA models,
               understand how the different modelling tools operate, and their limitations.
               Use the sidebar to simulate an ARIMA process. You can change the in-sample and the our-of-sample size. 
               Longer in-sample makes it easier for estimation and statistical tests to perform well. 
               Model parameters are generated randomly, unless set by the user. 
               All in-sample diagnostics are interactive. You can explore the ACF/PACH, set the differencing order for plots, 
               unit root tests, and information criteria. You get in-sample and out-of-sample accuracy (Root Mean Squared Error),
               with the option of using benchmarks of the true model with estimated parameters, 
               the theoretical model (known orders and parameters), and an automatically identified model
               using KPSS for the differencing order and AICc for the AR and MA orders (auto.arima).
               You can also consult in-sample residuals diagnostics, such as the Ljung-Box Q-test. 
               There is a clear separation between in-sample (fitting) and out-of-sample (testing) data,
               so that you can stress test the different tools and observe how they fail as the data generating
               model becomes more complex, or the in-sample statistics do not guard you against overfitting.")
           )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
      tags$h4("1. Data Setup"),
      helpText(tags$em("Select sample sizes, model orders, and click generate.")),
      wellPanel(
        sliderInput("n_train", "In-Sample size:", min = 20, max = 500, value = 60, step = 20),
        sliderInput("n_test", "Out-of-Sample size:", min = 20, max = 100, value = 20, step = 20),
        fluidRow(
          column(4, numericInput("p", "p", 1, min = 0, max = 10)),
          column(4, numericInput("d", "d", 0, min = 0, max = 3)),
          column(4, numericInput("q", "q", 1, min = 0, max = 10))
        ),
        checkboxInput("manual_mode", "Manual coefficients?", FALSE),
        conditionalPanel(
          condition = "input.manual_mode == true",
          helpText(tags$em("Separare multiple values with comma. Keep empty for 0 order.")),
          textInput("ar_manual", "AR Coeffs", "0.5"),
          textInput("ma_manual", "MA Coeffs", "0.4")
        ),
        actionButton("simulate", "Generate New Series", class = "btn-primary btn-block")
      ),
      hr(),
      tags$h4("2. Visual Exploration"),
      helpText(tags$em("Plots and statistics update with selected differencing.")),
      numericInput("diff_order", "Difference order for ACF/PACF:", 0, min = 0, max = 3),
      checkboxInput("plot_diff_series", "Plot differenced time series?", FALSE),
      checkboxInput("show_tests", "Show Unit Root Tests", FALSE),
      hr(),
      tags$h4("3. Your Model Identification"),
      helpText(tags$em("Select the model orders according to the specification strategy you prefer.")),
      fluidRow(
        column(4, numericInput("user_p", "p", 0, min = 0)),
        column(4, numericInput("user_d", "d", 0, min = 0)),
        column(4, numericInput("user_q", "q", 0, min = 0))
      ),
      checkboxInput("show_ic", "Show Information Criteria for User Model", FALSE),
      helpText(tags$em("Benchmark in- and out-of-sample accuracy.")),
      checkboxInput("show_truth", "Compare with True & Theoretical Models", FALSE),
      checkboxInput("show_auto", "Compare with Auto.Arima", FALSE),
      hr(),
      tags$h4("4. Forecast & Fit Visualization"),
      helpText(tags$em("Note: Forecast are only displayed when benchmark models are enabled in section 3.")),
      checkboxInput("show_fit", "Show In-Sample Fits", FALSE),
      checkboxInput("plot_user", "Plot Your Forecast", TRUE),
      checkboxInput("plot_true_est", "Plot True Order (Estimated Parameters)", FALSE),
      checkboxInput("plot_true_theo", "Plot Theoretical (Orders and Parameters Known)", FALSE),
      checkboxInput("plot_auto", "Plot Auto.Arima Forecast", FALSE),
      hr(),
      tags$h4("5. Residual Diagnostics"),
      checkboxInput("show_diag", "Run Ljung-Box & Residual Analysis", FALSE)
    ),
    
    mainPanel(width = 9,
      wellPanel(
        style = "background-color: #ffffff; border: 1px solid #ddd; text-align: center; ; padding: 0px 0px;",
        #tags$h5("Simulated Data Generating Process (DGP)"),
        uiOutput("arimaEq")
      ),
      plotOutput("tsPlot", height = "400px"),
      br(),
      fluidRow(
        column(6, plotOutput("acfPlot", height = "400px")),
        column(6, plotOutput("pacfPlot", height = "400px"))
      ),
      br(),
      uiOutput("dynamicResultsLayout"),
      hr(),
      conditionalPanel(
        condition = "input.show_diag == true",
        tags$h3("In-Sample Residual Diagnostics (Your Model)"),
        fluidRow(
          column(4, tags$h4("Ljung-Box Q-Test"), tableOutput("lbTable")),
          column(4, plotOutput("resAcf", height = "250px")),
          column(4, plotOutput("resPacf", height = "250px"))
        ),
        br(),
        fluidRow(
          column(6, plotOutput("resHist", height = "300px")),
          column(6, plotOutput("resQQ", height = "300px"))
        )
      )
    )
  ),
  fluidRow(
    column(12, hr(), div(style = "color: #888888; text-align: right; padding: 20px;", 
                         "Nikolaos Kourentzes | 2026-02-15"))
  )
)

server <- function(input, output, session) {
  
  observe({
    toggleState("plot_auto", input$show_auto)
    toggleState("plot_true_est", input$show_truth)
    toggleState("plot_true_theo", input$show_truth)
  })
  
  gen_stable_coefs <- function(order) {
    if (order == 0) return(numeric(0))
    coefs <- runif(order, -0.3, 0.3)
    if(sum(abs(coefs)) >= 1) coefs <- coefs / (sum(abs(coefs)) + 0.1)
    return(coefs)
  }
  
  sim_results <- eventReactive(input$simulate, {
    ar_vec <- if(input$manual_mode) as.numeric(unlist(strsplit(input$ar_manual, ","))) else gen_stable_coefs(input$p)
    ma_vec <- if(input$manual_mode) as.numeric(unlist(strsplit(input$ma_manual, ","))) else gen_stable_coefs(input$q)
    ts_data <- arima.sim(n = input$n_train + input$n_test, 
                         model = list(ar = ar_vec, ma = ma_vec, order = c(input$p, input$d, input$q)), sd = 1)
    list(train = subset(ts_data, end = input$n_train), test = subset(ts_data, start = input$n_train + 1),
         ar = ar_vec, ma = ma_vec, p = input$p, d = input$d, q = input$q, full = ts_data)
  }, ignoreNULL = FALSE)
  
  output$arimaEq <- renderUI({
    res <- sim_results()
    req(res)
    
    prefix <- "\\color{#377EB8FF}{AR}\\color{#999999FF}{I}\\color{#A65628FF}{MA}: "
    
    i_part <- if(res$d == 0) "" else paste0("\\color{#999999FF}{(1-B)^{", res$d, "}}")
    
    ar_terms <- if(length(res$ar) > 0) {
      paste0(sapply(seq_along(res$ar), function(i) {
        val <- res$ar[i]
        paste0(ifelse(val >= 0, " + ", " - "), "\\color{#377EB8FF}{\\mathbf{", abs(round(val, 3)), "}} y_{t-", i, "}")
      }), collapse = "")
    } else ""
    
    ma_terms <- if(length(res$ma) > 0) {
      paste0(sapply(seq_along(res$ma), function(i) {
        val <- res$ma[i]
        paste0(ifelse(val >= 0, " + ", " - "), "\\color{#A65628FF}{\\mathbf{", abs(round(val, 3)), "}} \\varepsilon_{t-", i, "}")
      }), collapse = "")
    } else ""
    
    lhs <- paste0(i_part, " y_t")
    
    # equation <- paste0(
    #   "$$\\begin{aligned} ", 
    #   prefix, " & ", lhs, " = ", ar_terms, " \\\\ ",
    #   " & ", ma_terms, " + \\varepsilon_t ",
    #   "\\end{aligned}$$"
    # )
    equation <- paste0(
      "$$", prefix, lhs, " = ", ar_terms, ma_terms, " + \\varepsilon_t", "$$"
    )
    
    equation <- gsub("=  \\+", "= ", equation)
    
    tags$div(
      style = "overflow-x: auto; overflow-y: hidden; width: 100%; padding: 10px 0;",
      withMathJax(helpText(equation))
    )
    # withMathJax(helpText(equation))
  })
  
  diag_data <- reactive({
    req(sim_results())
    x <- sim_results()$train
    if(input$diff_order > 0) x <- diff(x, differences = input$diff_order)
    x
  })
  
  models <- reactive({
    res <- sim_results()
    auto_mod <- if(input$show_auto) auto.arima(res$train) else NULL
    auto_label <- if(!is.null(auto_mod)) {
      ord <- arimaorder(auto_mod)
      paste0("Auto.Arima (", ord[1], ",", ord[2], ",", ord[3], ")")
    } else "Auto.Arima"
    
    list(
      user = tryCatch(Arima(res$train, order = c(input$user_p, input$user_d, input$user_q)), error = function(e) NULL),
      true_est = tryCatch(Arima(res$train, order = c(res$p, res$d, res$q), include.mean = FALSE), error = function(e) NULL),
      theo = tryCatch(Arima(res$train, order = c(res$p, res$d, res$q), fixed = c(res$ar, res$ma), include.mean = FALSE), error = function(e) NULL),
      auto = auto_mod,
      auto_label = auto_label
    )
  })
  
  render_custom_corr <- function(data, type="ACF") {
    corr_obj <- if(type=="ACF") acf(data, plot=FALSE) else pacf(data, plot=FALSE)
    df <- data.frame(lag = as.numeric(corr_obj$lag), val = as.numeric(corr_obj$acf))
    ci <- qnorm((1 + 0.95) / 2) / sqrt(length(data))
    ggplot(df, aes(x=lag, y=val)) +
      geom_segment(aes(xend=lag, yend=0), linewidth = 1.3, color="black") +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      geom_hline(yintercept = c(ci, -ci), color="blue", linetype="dashed", linewidth = 0.5) +
      scale_x_continuous(breaks = seq(0, max(df$lag), by=2)) +
      theme_minimal() + labs(title=paste("In-Sample", type, "| I =", input$diff_order), y=type) + ylim(-1,1)
  }
  
  output$acfPlot <- renderPlot({ render_custom_corr(diag_data(), "ACF") })
  output$pacfPlot <- renderPlot({ render_custom_corr(diag_data(), "PACF") })
  
  output$tsPlot <- renderPlot({
    res <- sim_results(); mods <- models()
    cols <- c("Observed (Future)"="black", "Your Model"="red", "True Order (Est)"="blue", "Theoretical"="green")
    cols[mods$auto_label] <- "purple"
    
    label_y <- max(res$full) + (0.15 * diff(range(res$full)))
    p <- autoplot(res$train, linewidth = 1.1) + 
      autolayer(res$test, series = "Observed (Future)", color = "black", linetype = "dashed", linewidth = 1.1) +
      geom_vline(xintercept = input$n_train, linetype="dotted", color="grey40") +
      annotate("text", x = input$n_train/2, y = label_y, label = "In-sample", fontface = "bold", size = 5) +
      annotate("text", x = input$n_train + (input$n_test/2), y = label_y, label = "Out-of-sample", fontface = "bold", size = 5) +
      theme_minimal() + labs(y = "Time series", x = "Time")
    
    if(input$plot_user && !is.null(mods$user)) {
      if(input$show_fit) p <- p + autolayer(fitted(mods$user), series = "Your Model")
      p <- p + autolayer(forecast(mods$user, h = input$n_test)$mean, series = "Your Model")
    }
    if(input$show_truth) {
      if(input$plot_true_est && !is.null(mods$true_est)) {
        if(input$show_fit) p <- p + autolayer(fitted(mods$true_est), series = "True Order (Est)")
        p <- p + autolayer(forecast(mods$true_est, h = input$n_test)$mean, series = "True Order (Est)")
      }
      if(input$plot_true_theo && !is.null(mods$theo)) {
        if(input$show_fit) p <- p + autolayer(fitted(mods$theo), series = "Theoretical")
        p <- p + autolayer(forecast(mods$theo, h = input$n_test)$mean, series = "Theoretical")
      }
    }
    if(input$show_auto && input$plot_auto && !is.null(mods$auto)) {
      if(input$show_fit) p <- p + autolayer(fitted(mods$auto), series = mods$auto_label)
      p <- p + autolayer(forecast(mods$auto, h = input$n_test)$mean, series = mods$auto_label)
    }
    p + scale_color_manual(values = cols)
  })
  
  output$dynamicResultsLayout <- renderUI({
    fluidRow(
      column(4, 
             if(input$show_ic) tagList(
               h4("Information Criteria for Your Model"), 
               tableOutput("icTable"),
               helpText(tags$em("Note: Compare models only with the same differencing order (d)."))
             ) else NULL,
             h4("Forecast Accuracy (RMSE)"), 
             tableOutput("rmseTable"),
             helpText(tags$em("Best (★) is based on the lowest Out-of-Sample RMSE. True & Theoretical models are unknown in practice."))
      ),
      column(4, 
             conditionalPanel(condition = "input.show_tests == true", 
                              h4(paste0("In-Sample Unit Root Tests (I = ", input$diff_order, ")")), 
                              tableOutput("statTestsTable"),
                              helpText(tags$b("Interpretation:")),
                              helpText("1. For DF, ADF, and PP: Stationarity if Stat < Critical Value (Null = Unit Root)."),
                              helpText("2. For KPSS: Stationarity if Stat < Critical Value (Null = Stationarity).")
             )
      ),
      column(4, 
             if(input$plot_diff_series) tagList(
               h4(paste0("In-Sample Series (I = ", input$diff_order, ")")), 
               plotOutput("diffTsPlot", height = "250px")
             ) else NULL)
    )
  })
  
  output$icTable <- renderTable({
    m <- models()$user
    if(is.null(m)) return(NULL)
    data.frame(Criterion = c("AIC", "AICc", "BIC"), Value = c(m$aic, m$aicc, m$bic))
  }, digits = 2)
  
  output$rmseTable <- renderTable({
    res <- sim_results(); mods <- models(); h = input$n_test; results <- list()
    calc_metrics <- function(m, label) {
      if(is.null(m)) return(NULL)
      fc <- forecast(m, h=h); rmse_out <- sqrt(mean((res$test - fc$mean)^2)); rmse_in <- sqrt(mean(residuals(m)^2))
      data.frame(Model=label, `In-Sample`=rmse_in, `Out-of-Sample`=rmse_out)
    }
    results[[1]] <- calc_metrics(mods$user, "Your Model")
    if(input$show_truth) { results[[2]] <- calc_metrics(mods$true_est, "True (Est)"); results[[3]] <- calc_metrics(mods$theo, "Theoretical") }
    if(input$show_auto) results[[4]] <- calc_metrics(mods$auto, mods$auto_label)
    df <- do.call(rbind, results)
    df$Best <- ifelse(df$Out.of.Sample == min(df$Out.of.Sample, na.rm=TRUE), "★", "")
    df
  }, digits = 3)
  
  output$statTestsTable <- renderTable({
    x <- diag_data()
    t_df_n <- ur.df(x, type="none", lags=0); t_df_d <- ur.df(x, type="drift", lags=0); t_df_t <- ur.df(x, type="trend", lags=0)
    t_adf_n <- ur.df(x, type="none", lags=1); t_adf_d <- ur.df(x, type="drift", lags=1); t_adf_t <- ur.df(x, type="trend", lags=1)
    t_pp <- ur.pp(x, type="Z-tau", model="trend"); t_kpss <- ur.kpss(x, type="tau")
    get_row <- function(obj, name) {
      stat <- obj@teststat[1]; crit <- obj@cval[1,2] 
      data.frame(Test=name, Statistic=round(stat,3), `Crit Value (5%)`=round(crit,3), Outcome=if(stat < crit) "Stationary" else "Non-stationary")
    }
    rbind(get_row(t_df_n, "DF (None)"), get_row(t_df_d, "DF (Drift)"), get_row(t_df_t, "DF (Trend)"),
          get_row(t_adf_n, "ADF (None)"), get_row(t_adf_d, "ADF (Drift)"), get_row(t_adf_t, "ADF (Trend)"),
          get_row(t_pp, "PP Test"), get_row(t_kpss, "KPSS"))
  })
  
  output$diffTsPlot <- renderPlot({ 
    autoplot(diag_data(), size = 1) + theme_minimal() + 
      labs(y = "Value", title = "") 
  })
  
  output$lbTable <- renderTable({
    m <- models()$user; if(is.null(m)) return(NULL); res <- residuals(m)
    bt10 <- Box.test(res, lag=10, type="Ljung-Box"); bt20 <- Box.test(res, lag=20, type="Ljung-Box")
    data.frame(Lag=c(10, 20), Statistic=round(c(bt10$statistic, bt20$statistic), 3), `p-value`=round(c(bt10$p.value, bt20$p.value), 3), Outcome=ifelse(c(bt10$p.value, bt20$p.value) > 0.05, "White Noise", "Serial Correlation"))
  })
  output$resAcf <- renderPlot({ if(!is.null(models()$user)) render_custom_corr(residuals(models()$user), "ACF") })
  output$resPacf <- renderPlot({ if(!is.null(models()$user)) render_custom_corr(residuals(models()$user), "PACF") })
  output$resHist <- renderPlot({ if(is.null(models()$user)) return(NULL); res <- as.numeric(residuals(models()$user)); ggplot(data.frame(res), aes(x=res)) + geom_histogram(aes(y=..density..), bins=15, fill="steelblue", alpha=0.7) + stat_function(fun=dnorm, args=list(mean=mean(res), sd=sd(res)), color="red", size=1) + theme_minimal() + ggtitle("In-Sample Hist. of Residuals") })
  output$resQQ <- renderPlot({ if(is.null(models()$user)) return(NULL); ggplot(data.frame(res = as.numeric(residuals(models()$user))), aes(sample=res)) + stat_qq() + stat_qq_line(color="red") + theme_minimal() + ggtitle("In-Sample Normal Q-Q Plot") })
}

shinyApp(ui, server)