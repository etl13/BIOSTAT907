# ========================================================================
# SHINY APP — Two-Stage UMVUE + MLE With CI + PMF + P-values
# Fully updated version with:
#  • estimator-ordering p-values (Jung, Lecture 3)
#  • tabs for estimates, CIs, p-values, plots
#  • PMF hide/show toggle
# ========================================================================

library(shiny)
library(tidyverse)
library(DT)
library(bslib)

# ========================================================
# Helper function: UMVUE, MLE, and f(m,s|p) table
# ========================================================
umvue_mle_table_with_p <- function(n1, n2, a1, a = NULL,
                                   p_vec = c(0.1,0.2,0.3,0.4,0.5)) {
  N <- n1 + n2
  
  results <- data.frame(
    m = integer(),
    s = integer(),
    UMVUE = numeric(),
    MLE   = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add columns for f(m,s|p)
  for (p in p_vec) {
    results[[paste0("f_p", gsub("\\.", "", as.character(p)))]] <- numeric()
  }
  
  # ---------- Early stop rows ----------
  for (x1 in 0:a1) {
    umv <- x1 / n1
    mle <- x1 / n1
    newrow <- data.frame(m=1, s=x1, UMVUE=umv, MLE=mle)
    
    # PMFs
    for (p in p_vec) {
      prob <- dbinom(x1, n1, p)
      newrow[[paste0("f_p", gsub("\\.", "", as.character(p)))]] <- prob
    }
    
    results <- rbind(results, newrow)
  }
  
  # ---------- Continuation rows ----------
  for (t in (a1+1):N) {
    x1_min <- max(0, t - n2)
    x1_max <- min(n1, t)
    x1_cond_min <- max(a1+1, x1_min)
    
    if (x1_cond_min <= x1_max) {
      x1_vals <- x1_cond_min:x1_max
      
      # UMVUE expectation (Rao-Blackwell)
      logp <- lchoose(n1, x1_vals) + 
        lchoose(n2, t - x1_vals) -
        lchoose(N, t)
      logp <- logp - max(logp)
      probs <- exp(logp)
      probs <- probs / sum(probs)
      E_x1 <- sum(x1_vals * probs)
      umv <- E_x1 / n1
      mle <- t / N
      
      newrow <- data.frame(m=2, s=t, UMVUE=umv, MLE=mle)
      
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
  
  # Normalize all f(m,s|p)
  for (p in p_vec) {
    col <- paste0("f_p", gsub("\\.", "", as.character(p)))
    results[[col]] <- results[[col]] / sum(results[[col]])
  }
  
  rownames(results) <- NULL
  return(results)
}

# ========================================================
# Exact CI (your existing function, unchanged)
# ========================================================
ci_two_stage_exact <- function(n1, n2, a1,
                               obs_m, obs_s,
                               estimator=c("MLE","UMVUE"),
                               alpha=0.05,
                               tol=1e-8,
                               maxiter=100) {
  estimator <- match.arg(estimator)
  rows <- umvue_mle_table_with_p(n1, n2, a1)
  idx_obs <- which(rows$m == obs_m & rows$s == obs_s)
  if (length(idx_obs) != 1) stop("Observed (m,s) not found.")
  
  obs_val <- if (estimator=="UMVUE") rows$UMVUE[idx_obs] else rows$MLE[idx_obs]
  
  row_probs_given_p <- function(p) {
    probs <- numeric(nrow(rows))
    for (j in seq_len(nrow(rows))) {
      if (rows$m[j] == 1) {
        x1 <- rows$s[j]
        probs[j] <- dbinom(x1, n1, p)
      } else {
        t <- rows$s[j]
        x1_min <- max(a1+1, t - n2)
        x1_max <- min(n1, t)
        if (x1_min <= x1_max) {
          x1_vals <- x1_min:x1_max
          probs[j] <- sum(dbinom(x1_vals,n1,p)*dbinom(t-x1_vals,n2,p))
        }
      }
    }
    s <- sum(probs)
    if (s <= 0) return(probs * 0)
    probs / s
  }
  
  cum_le <- function(p) {
    probs <- row_probs_given_p(p)
    est_vals <- if (estimator=="UMVUE") rows$UMVUE else rows$MLE
    sum(probs[est_vals <= obs_val+1e-12])
  }
  cum_ge <- function(p) {
    probs <- row_probs_given_p(p)
    est_vals <- if (estimator=="UMVUE") rows$UMVUE else rows$MLE
    sum(probs[est_vals >= obs_val-1e-12])
  }
  
  # lower CI
  f_low <- function(p) cum_ge(p) - alpha/2
  p_lower <- NA_real_
  if (f_low(0)>=0) p_lower <- 0 else if (f_low(1)>=0) {
    p_lower <- uniroot(f_low,lower=0,upper=1,tol=tol)$root
  }
  
  # upper CI
  g_up <- function(p) cum_le(p) - alpha/2
  p_upper <- NA_real_
  if (g_up(1)>=0) p_upper <- 1 else if (g_up(0)>=0) {
    p_upper <- uniroot(g_up,lower=0,upper=1,tol=tol)$root
  }
  
  c(lower=p_lower, upper=p_upper)
}

# ========================================================
# NEW: Correct estimator-ordering p-value (Jung lecture)
# ========================================================
p_value_two_stage <- function(n1, n2, a1, m_obs, s_obs, p0,
                              estimator=c("UMVUE","MLE")) {
  estimator <- match.arg(estimator)
  
  rows <- umvue_mle_table_with_p(n1, n2, a1, p_vec=p0)
  pmf_col <- paste0("f_p", gsub("\\.","",as.character(p0)))
  
  idx_obs <- which(rows$m == m_obs & rows$s == s_obs)
  if (length(idx_obs) != 1) stop("Observed (m,s) not found.")
  
  est_vals <- if (estimator=="UMVUE") rows$UMVUE else rows$MLE
  fvec <- rows[[pmf_col]]
  fvec <- fvec / sum(fvec)
  
  obs_est <- est_vals[idx_obs]
  
  # Upper-tail p-value: P(est >= obs_est)
  sum(fvec[est_vals >= (obs_est - 1e-12)])
}

# ========================================================
# Combine everything into complete table
# ========================================================
two_stage_ci_table <- function(n1, n2, a1, a, alpha=0.05,
                               p_values=c(0.1,0.2,0.3,0.4,0.5),
                               p0=0.2) {
  
  out <- umvue_mle_table_with_p(n1, n2, a1, a, p_values)
  
  out$MLE_lower <- NA_real_
  out$MLE_upper <- NA_real_
  out$UMVUE_lower <- NA_real_
  out$UMVUE_upper <- NA_real_
  out$MLE_pvalue <- NA_real_
  out$UMVUE_pvalue <- NA_real_
  
  for (i in seq_len(nrow(out))) {
    m_i <- out$m[i]
    s_i <- out$s[i]
    
    ci_mle <- ci_two_stage_exact(n1, n2, a1, m_i, s_i, "MLE", alpha)
    ci_umvue <- ci_two_stage_exact(n1, n2, a1, m_i, s_i, "UMVUE", alpha)
    
    out$MLE_lower[i] <- ci_mle[1]
    out$MLE_upper[i] <- ci_mle[2]
    out$UMVUE_lower[i] <- ci_umvue[1]
    out$UMVUE_upper[i] <- ci_umvue[2]
    
    out$MLE_pvalue[i] <- p_value_two_stage(n1, n2, a1, m_i, s_i, p0, "MLE")
    out$UMVUE_pvalue[i] <- p_value_two_stage(n1, n2, a1, m_i, s_i, p0, "UMVUE")
  }
  
  out
}

# ========================================================
# PMF Plot helpers
# ========================================================
plot_pmf1 <- function(check) {
  check %>%
    pivot_longer(cols=starts_with("f"), names_to="func", values_to="val") %>%
    ggplot(aes(UMVUE, val)) +
    geom_col() +
    labs(x="", y="UMVUE", title="PMF (UMVUE ordering)") +
    theme_bw() +
    facet_wrap(~func)
}

plot_pmf2 <- function(check) {
  check %>%
    pivot_longer(cols=starts_with("f"), names_to="func", values_to="val") %>%
    ggplot(aes(MLE, val)) +
    geom_col() +
    labs(x="", y="MLE", title="PMF (MLE ordering)") +
    theme_bw() +
    facet_wrap(~func)
}

# ========================================================
# SERVER
# ========================================================
server <- function(input, output, session) {
  
  table_data <- eventReactive(input$run, {
    n1 <- input$n1; n2 <- input$n2
    a1 <- input$a1; a  <- input$a
    alpha <- input$alpha
    p_values <- as.numeric(strsplit(input$pvals, ",")[[1]])
    p0 <- input$p0
    
    two_stage_ci_table(n1, n2, a1, a, alpha, p_values, p0)
  })
  
  # ------------------ ESTIMATES TABLE -----------------
  output$estimates_table <- renderDT({
    req(table_data())
    df <- table_data()
    pmf_cols <- grep("^f_", names(df), value=TRUE)
    base_cols <- c("m","s","UMVUE","MLE")
    
    if (isTRUE(input$show_pmf))
      df2 <- df[, c(base_cols, pmf_cols)]
    else
      df2 <- df[, base_cols]
    
    datatable(df2, options=list(pageLength=20, scrollX=TRUE))
  })
  
  # ------------------ CI TABLE -------------------------
  output$ci_table <- renderDT({
    req(table_data())
    df <- table_data()
    pmf_cols <- grep("^f_", names(df), value=TRUE)
    base_cols <- c("m","s","UMVUE_lower","UMVUE_upper",
                   "MLE_lower","MLE_upper")
    
    if (isTRUE(input$show_pmf))
      df2 <- df[, c(base_cols, pmf_cols)]
    else
      df2 <- df[, base_cols]
    
    datatable(df2, options=list(pageLength=20, scrollX=TRUE))
  })
  
  # ------------------ P-VALUE TABLE ---------------------
  output$pvalue_table <- renderDT({
    req(table_data())
    df <- table_data()
    
    pmf_cols <- grep("^f_", names(df), value=TRUE)
    base_cols <- c("m","s","UMVUE","MLE",
                   "UMVUE_pvalue","MLE_pvalue")
    
    if (isTRUE(input$show_pmf))
      df2 <- df[, c(base_cols, pmf_cols)]
    else
      df2 <- df[, base_cols]
    
    datatable(df2, options=list(pageLength=20, scrollX=TRUE))
  })
  
  # ------------------ PMF PLOTS -------------------------
  output$results_plot1 <- renderPlot({
    req(table_data())
    plot_pmf1(table_data())
  })
  
  output$results_plot2 <- renderPlot({
    req(table_data())
    plot_pmf2(table_data())
  })
  
  # ------------------ DOWNLOAD --------------------------
  output$downloadData <- downloadHandler(
    filename = function() paste0("two_stage_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(table_data(), file, row.names=FALSE)
  )
}

