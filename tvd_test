# Input must be a data frame with ea. row corresponding to a CpG that includes columns specifying direction (5' or 3'), the bin that it belongs to and MP value (%).
# Default column header for the required columns are: "Dir", "Bin" and "MP"
# Bin size and bound can be set.

library(data.table)
library(dplyr)

bin_size=10; bound=2000

tvd_test <- function(df, value_col = "MP", dir_col = "Dir", bin_col = "Bin", n_perm = 1000) {
  calc_tvd <- function(x, y, breaks = NULL) {
    if (is.null(breaks)) {
      min_val <- min(c(x, y))
      max_val <- max(c(x, y))
      breaks <- seq(min_val, max_val, length.out = 30)
    }
    px <- hist(x, breaks = breaks, plot = FALSE)$counts
    py <- hist(y, breaks = breaks, plot = FALSE)$counts
    px <- px / sum(px)
    py <- py / sum(py)
    0.5 * sum(abs(px - py))
  }
  df[[dir_col]] <- as.character(df[[dir_col]])
  unique_bins <- unique(df[[bin_col]])
  results_list <- lapply(unique_bins, function(bin_val) {
    subdf <- subset(df, df[[bin_col]] == bin_val)
    subdf[[dir_col]] <- as.character(subdf[[dir_col]])

    group_vals <- unique(subdf[[dir_col]])
    if(length(group_vals) != 2) {
      cat(sprintf("Bin: %s - Skipped (does not have exactly two groups)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }

    x1 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[1]]
    x2 <- subdf[[value_col]][subdf[[dir_col]] == group_vals[2]]

    if(length(x1) == 0 || length(x2) == 0) {
      cat(sprintf("Bin: %s - Skipped (empty group)\n", bin_val))
      return(data.frame(Bin = bin_val, TVD = NA, p_value = NA, Larger_Dir = NA))
    }

    obs_tvd <- calc_tvd(x1, x2)
    pooled <- c(x1, x2)
    n1 <- length(x1)

    perm_tvd <- replicate(n_perm, {
      perm_labels <- sample(rep(group_vals, c(n1, length(x2))))
      calc_tvd(pooled[perm_labels == group_vals[1]], pooled[perm_labels == group_vals[2]])
    })

    p_value <- mean(perm_tvd >= obs_tvd)

    mean1 <- mean(x1)
    mean2 <- mean(x2)

    larger_dir <- ifelse(mean1 > mean2, group_vals[1], group_vals[2])

    cat(sprintf("Bin: %s, TVD: %.4f, p-value: %.4f, Larger Dir: %s\n", bin_val, obs_tvd, p_value, larger_dir))

    data.frame(Bin = bin_val, TVD = obs_tvd, p_value = p_value, Larger_Dir = larger_dir)
  })
  do.call(rbind, results_list)
}
