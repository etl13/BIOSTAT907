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

# p-value based on estimator-ordering per Jung (lecture ch.3)
# n1, n2, a1 : design
# m_obs, s_obs : observed stopping stage and cumulative successes
# p0 : null value (H0: p = p0)
# estimator : "UMVUE" or "MLE"
p_value_two_stage <- function(n1, n2, a1, m_obs, s_obs, p0,
                              estimator = c("UMVUE", "MLE")) {
  estimator <- match.arg(estimator)
  # get table with PMF column for p0
  # umvue_mle_table_with_p normalizes f(...) column so it sums to 1
  rows <- umvue_mle_table_with_p(n1 = n1, n2 = n2, a1 = a1,
                                 a = NULL,
                                 p_vec = p0)
  # name of pmf column produced by umvue_mle_table_with_p
  pmf_col <- paste0("f_p", gsub("\\.", "", as.character(p0)))
  if (!pmf_col %in% names(rows)) stop("PMF column not found in rows.")
  
  # find observed row index
  idx_obs <- which(rows$m == m_obs & rows$s == s_obs)
  if (length(idx_obs) != 1) stop("Observed (m,s) not found in table rows.")
  
  # pick estimator values and the f vector
  est_vals <- if (estimator == "UMVUE") rows$UMVUE else rows$MLE
  fvec <- rows[[pmf_col]]
  # numerical safety: ensure fvec sums to 1
  if (sum(fvec) <= 0) stop("Computed PMF sums to 0 at p0.")
  fvec <- fvec / sum(fvec)
  
  # observed estimator
  obs_est <- est_vals[idx_obs]
  
  # p-value (upper-tail): sum f(i,j|p0) for estimator >= obs_est
  # include equality (>=) per definition in slides
  pval <- sum(fvec[est_vals >= (obs_est - .Machine$double.eps)])
  
  return(as.numeric(pval))
}



# --- Combine all
two_stage_ci_table <- function(n1, n2, a1, a, alpha = 0.05, p_values = c(0.1,0.2,0.3,0.4,0.5)) {
  out <- umvue_mle_table_with_p(n1, n2, a1, alpha, p_values)
  out$MLE_lower <- NA_real_
  out$MLE_upper <- NA_real_
  out$UMVUE_lower <- NA_real_
  out$UMVUE_upper <- NA_real_
  out$MLE_pvalue <- NA_real_
  out$UMVUE_pvalue <- NA_real_
  

  for (i in seq_len(nrow(out))) {
    ci_mle <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "MLE", alpha)
    ci_umvue <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "UMVUE", alpha)
    out$MLE_lower[i] <- ci_mle[1]
    out$MLE_upper[i] <- ci_mle[2]
    out$UMVUE_lower[i] <- ci_umvue[1]
    out$UMVUE_upper[i] <- ci_umvue[2]
  }
  for (i in seq_len(nrow(out))) {
    ci_mle <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "MLE", alpha)
    ci_umvue <- ci_two_stage_exact(n1, n2, a1, out$m[i], out$s[i], "UMVUE", alpha)
    out$MLE_lower[i] <- ci_mle[1]
    out$MLE_upper[i] <- ci_mle[2]
    out$UMVUE_lower[i] <- ci_umvue[1]
    out$UMVUE_upper[i] <- ci_umvue[2]
    
    # NEW: exact p-values using Jung-style ordering (use p0 from input)
    # Here you must decide which p0 to use for the test; typically that's an input value.
    # Example: assume you add 'p0' as an input or reuse the first element of p_values.
    # I'll use 'p0' variable below — make sure it exists in your scope.
    out$MLE_pvalue[i] <- p_value_two_stage(n1, n2, a1, out$m[i], out$s[i], p0, "MLE")
    out$UMVUE_pvalue[i] <- p_value_two_stage(n1, n2, a1, out$m[i], out$s[i], p0, "UMVUE")
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
    a  <- input$a
    p0<-input$p0
    alpha <- input$alpha
    p_values <- as.numeric(strsplit(input$pvals, ",")[[1]])
    
    two_stage_ci_table(n1, n2, a1, a, alpha, p_values)
  })
  
  # ----------------------------------------------------
  # ESTIMATE TABLE: m, s, UMVUE, MLE, and PMF columns
  # ----------------------------------------------------
  output$estimates_table <- renderDT({
    req(table_data())
    df <- table_data()
    pmf_cols <- grep("^f_", names(df), value = TRUE)
    est_df <- df[, c("m", "s", "UMVUE", "MLE", pmf_cols)]
    
    datatable(est_df, options = list(pageLength = 20, scrollX = TRUE))
  })
  
  # ----------------------------------------------------
  # CI TABLE: m, s, UMVUE CIs, MLE CIs, and PMF columns
  # ----------------------------------------------------
  output$ci_table <- renderDT({
    req(table_data())
    df <- table_data()
    pmf_cols <- grep("^f_", names(df), value = TRUE)
    
    ci_df <- df[, c("m", "s",
                    "UMVUE_lower", "UMVUE_upper",
                    "MLE_lower", "MLE_upper",
                    pmf_cols)]
    
    datatable(ci_df, options = list(pageLength = 20, scrollX = TRUE))
  })
  
  output$pvalue_table <- renderDT({
    req(table_data())
    df <- table_data()
    pmf_cols <- grep("^f_", names(df), value = TRUE)
    
    base_cols <- c("m", "s", "UMVUE", "MLE", "UMVUE_pvalue", "MLE_pvalue")
    if (isTRUE(input$show_pmf)) {
      pv_df <- df[, c(base_cols, pmf_cols)]
    } else {
      pv_df <- df[, base_cols]
    }
    datatable(pv_df, options = list(pageLength = 20, scrollX = TRUE))
  })
  
  
  # ---- PMF plots ----
  output$results_plot1 <- renderPlot({
    req(table_data())
    plot_pmf1(table_data())
  })
  
  output$results_plot2 <- renderPlot({
    req(table_data())
    plot_pmf2(table_data())
  })
  
  # ---- Download ----
  output$downloadData <- downloadHandler(
    filename = function() paste0("two_stage_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(table_data(), file, row.names = FALSE)
  )
}

