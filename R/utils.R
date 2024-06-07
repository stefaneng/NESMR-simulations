# Extracts lower triangular into FROM-TO format
extract_lower_triangular <- function(mat, prefix = "") {
  if (is.null(colnames(mat))) {
    col_names <- paste0("V", seq_len(ncol(mat)))
  } else {
    col_names <- colnames(mat)
  }

  low_tri_idx <- which(lower.tri(mat), arr.ind = TRUE)
  low_tri_arrow <- paste0(prefix, col_names[low_tri_idx[, 1]], '_', col_names[low_tri_idx[, 2]])
  setNames(mat[lower.tri(mat)], low_tri_arrow)
}


