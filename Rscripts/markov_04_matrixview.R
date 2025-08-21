
matrix_view <- function(mat, as_percent = FALSE, digits = 2) {
  # Create base matrix
  base <- mat
  
  # Optionally convert to percentage
  if (as_percent) {
    base <- round(base * 100, digits)
  }
  
  # Add row totals
  base_with_rowsum <- cbind(base, RowTotal = rowSums(base))
  
  # Add column totals (including the new row total column)
  full <- rbind(base_with_rowsum, ColTotal = colSums(base_with_rowsum))
  
  return(full)
}

matrix_view(transition_matrices[["FT24_1_1"]])

matrix_view(transition_matrices[["TR23_1_1"]])

matrix_view(transition_probabilities[["TR20_4_5"]])

matrix_view(transition_probabilities[["FT20_7_8"]])

