# =============================================================
# Shiny App: Two-Stage UMVUE and MLE with f(m,s|p) and CI Table
# =============================================================

library(shiny)
library(DT)
library(bslib)
library(tidyverse)

# -------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------
## -----------------------------------------------------------
## Version 5: UMVUE, MLE + normalized f(m,s|p)
## -----------------------------------------------------------

umvue_mle_table_with_p <- function(n1, n2, a1, a = NULL, 
                                   p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  N <- n1 + n2
  results <- data.frame(
    m = integer(),
    s = integer(),
    UMVUE = numeric(),
    MLE = numeric(),
    stringsAsFactors = FALSE
  )
  
  ## probability columns
  for (p in p_vec) {
    results[[paste0("f_p", gsub("\\.", "", as.character(p)))]] <- numeric()
  }
  
  ## 1) Early-stop rows
  for (x1 in 0:a1) {
    umv <- x1 / n1
    mle <- x1 / n1
    newrow <- data.frame(m = 1, s = x1, UMVUE = umv, MLE = mle)
    
    for (p in p_vec) {
      prob <- dbinom(x1, n1, p)
      newrow[[paste0("f_p", gsub("\\.", "", as.character(p)))]] <- prob
    }
    
    results <- rbind(results, newrow)
  }
  
  ## 2) Continuation rows
  for (t in (a1 + 1):N) {
    x1_min <- max(0, t - n2)
    x1_max <- min(n1, t)
    x1_cond_min <- max(x1_min, a1 + 1)
    
    if (x1_cond_min <= x1_max) {
      x1_vals <- seq(x1_cond_min, x1_max)
      
      # expectation for UMVUE (Rao-Blackwellization)
      logp <- lchoose(n1, x1_vals) + lchoose(n2, t - x1_vals) - lchoose(N, t)
      logp <- logp - max(logp)
      probs <- exp(logp)
      probs <- probs / sum(probs)
      E_x1 <- sum(x1_vals * probs)
      umv <- E_x1 / n1
      mle <- t / N
      
      newrow <- data.frame(m = 2, s = t, UMVUE = umv, MLE = mle)
      
      for (p in p_vec) {
        prob <- 0
        for (x1 in x1_cond_min:x1_max) {
          x2 <- t - x1
          prob <- prob + dbinom(x1, n1, p) * dbinom(x2, n2, p)
        }
        newrow[[paste0("f_p", gsub("\\.", "", as.character(p)))]] <- prob
      }
      
      results <- rbind(results, newrow)
    }
  }
  
  ## normalize f(m,s|p) so columns sum to 1
  for (p in p_vec) {
    col <- paste0("f_p", gsub("\\.", "", as.character(p)))
    results[[col]] <- results[[col]] / sum(results[[col]])
  }
  
  rownames(results) <- NULL
  return(results)
}
ci_two_stage_exact <- function(n1, n2, a1,
                               obs_m, obs_s,
                               estimator = c("MLE", "UMVUE"),
                               alpha = 0.05,
                               tol = 1e-8,
                               maxiter = 100) {
  estimator <- match.arg(estimator)
  N <- n1 + n2
  
  # build the compact rows once
  rows <- umvue_mle_table_with_p(n1, n2, a1, a)
  
  # find the observed row
  idx_obs <- which(rows$m == obs_m & rows$s == obs_s)
  if (length(idx_obs) != 1) stop("Observed (m,s) not found in table rows.")
  
  # observed estimator value
  obs_val <- if (estimator == "UMVUE") rows$UMVUE[idx_obs] else rows$MLE[idx_obs]
  
  # function to compute probability vector of rows under a given p
  row_probs_given_p <- function(p) {
    probs <- numeric(nrow(rows))
    for (j in seq_len(nrow(rows))) {
      if (rows$m[j] == 1L) {
        x1 <- rows$s[j]
        probs[j] <- dbinom(x1, n1, p)
      } else {
        t <- rows$s[j]
        x1_min <- max(a1 + 1, t - n2)
        x1_max <- min(n1, t)
        if (x1_min > x1_max) {
          probs[j] <- 0
        } else {
          x1_vals <- x1_min:x1_max
          probs_j <- dbinom(x1_vals, n1, p) * dbinom(t - x1_vals, n2, p)
          probs[j] <- sum(probs_j)
        }
      }
    }
    # normalize to sum to 1 (protect against tiny numeric error)
    s <- sum(probs)
    if (s <= 0) return(probs * 0)
    probs / s
  }
  
  # cumulative functions (<= obs_val, >= obs_val)
  cum_le <- function(p) {
    probs <- row_probs_given_p(p)
    est_vals <- if (estimator == "UMVUE") rows$UMVUE else rows$MLE
    sum(probs[est_vals <= obs_val + .Machine$double.eps])  # include equality
  }
  cum_ge <- function(p) {
    probs <- row_probs_given_p(p)
    est_vals <- if (estimator == "UMVUE") rows$UMVUE else rows$MLE
    sum(probs[est_vals >= obs_val - .Machine$double.eps])
  }
  
  ## find lower bound p_L solving cum_ge(p) = alpha/2 (smallest p with cum_ge >= alpha/2)
  f_low <- function(p) cum_ge(p) - (alpha / 2)
  f0_low <- f_low(0.0)
  f1_low <- f_low(1.0)
  p_lower <- NA_real_
  if (f0_low >= 0) {
    p_lower <- 0
  } else if (f1_low < 0) {
    # cannot reach alpha/2 even at p=1
    p_lower <- NA_real_
  } else {
    # bracketed: use uniroot
    try({
      root <- uniroot(f_low, lower = 0, upper = 1,
                      tol = tol, maxiter = maxiter)
      p_lower <- root$root
    }, silent = TRUE)
    if (is.na(p_lower)) {
      # fallback grid
      pseq <- seq(0, 1, length.out = 2001)
      vals <- sapply(pseq, function(pp) cum_ge(pp))
      idx <- which(vals >= alpha / 2)
      p_lower <- if (length(idx) == 0) NA_real_ else pseq[min(idx)]
    }
  }
  
  ## find upper bound p_U solving cum_le(p) = alpha/2 (largest p with cum_le >= alpha/2)
  g_up <- function(p) cum_le(p) - (alpha / 2)
  g0_up <- g_up(0.0)
  g1_up <- g_up(1.0)
  p_upper <- NA_real_
  if (g0_up < 0) {
    # even at p=0 lower-tail below alpha/2 -> no valid upper
    p_upper <- NA_real_
  } else if (g1_up >= 0) {
    p_upper <- 1
  } else {
    # bracket: g(0)>=0 and g(1)<0 so root exists
    try({
      root2 <- uniroot(g_up, lower = 0, upper = 1,
                       tol = tol, maxiter = maxiter)
      p_upper <- root2$root
    }, silent = TRUE)
    if (is.na(p_upper)) {
      # fallback grid
      pseq <- seq(0, 1, length.out = 2001)
      vals <- sapply(pseq, function(pp) cum_le(pp))
      idx <- which(vals >= alpha / 2)
      p_upper <- if (length(idx) == 0) NA_real_ else pseq[max(idx)]
    }
  }
  
  return(c(lower = p_lower, upper = p_upper))
}



# --- Combine all
two_stage_ci_table <- function(n1, n2, a1, a, alpha = 0.05, p_values = c(0.1,0.2,0.3,0.4,0.5)) {
  out <- umvue_mle_table_with_p(n1, n2, a1, alpha, p_values)
  out$MLE_lower <- NA_real_
  out$MLE_upper <- NA_real_
  out$UMVUE_lower <- NA_real_
  out$UMVUE_upper <- NA_real_
  
  for (i in seq_len(nrow(out))) {
    ci_mle <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "MLE", alpha)
    ci_umvue <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "UMVUE", alpha)
    out$MLE_lower[i] <- ci_mle[1]
    out$MLE_upper[i] <- ci_mle[2]
    out$UMVUE_lower[i] <- ci_umvue[1]
    out$UMVUE_upper[i] <- ci_umvue[2]
  }
  
  out
}


plot_pmf1<-function(check){
  
  check2<-check%>%
    pivot_longer(cols = c(starts_with("f")), names_to = "func", values_to = "val")
  
  ggplot(check2, aes(UMVUE, val))+
    geom_col()+
    labs(x = "", y = "UMVUE")+
    theme_bw()+
    facet_wrap(~func)
}
plot_pmf2<-function(check){
  check2<-check%>%
    pivot_longer(cols = c(starts_with("f")), names_to = "func", values_to = "val")
  
  ggplot(check2, aes(MLE, val))+
    geom_col()+
    labs(x = "", y = "MLE")+
    theme_bw()+
    facet_wrap(~func)
}




# =============================================================
# Server
# =============================================================
server <- function(input, output, session) {
  
  # ---- Reactive table computation ----
  table_data <- eventReactive(input$run, {
    n1 <- input$n1
    n2 <- input$n2
    a1 <- input$a1
    a <- input$a
    alpha <- input$alpha
    p_values <- as.numeric(strsplit(input$pvals, ",")[[1]])
    
    tbl <- two_stage_ci_table(n1, n2, a1, a, alpha, p_values)
    tbl
  })
  
  # ---- Render table ----
  output$results_table <- renderDT({
    req(table_data())
    datatable(table_data(), options = list(pageLength = 15, scrollX = TRUE))
  })
  
  # ---- Render plot ----
  output$results_plot1 <- renderPlot({
    req(table_data())
    plot_pmf1(table_data())
  })
  
  # ---- Render plot ----
  output$results_plot2<- renderPlot({
    req(table_data())
    plot_pmf2(table_data())
  })
  
  # ---- Download handler ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("two_stage_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(table_data(), file, row.names = FALSE)
    }
  )
}


