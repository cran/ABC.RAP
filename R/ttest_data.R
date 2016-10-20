#' applying t-test analysis
#' 
#' @description 
#' This function applies "two.sided", unequal variance Student's t-test analysis for each probe comparing cases and controls. A cutoff for p-values can be entered to minimise multiple testing bias to filter insignificant p-values.
#' 
#' @param x The filtered 450k probes from filter_data() function
#' @param cases_column_1 The first column (column number) for cases in the filtered dataset 
#' @param cases_column_n The last column (column number) for cases in the filtered dataset
#' @param controls_column_1 The first column (column number) for controls in the filtered dataset
#' @param controls_column_n The last column (column number) for controls in the filtered dataset
#' @param ttest_cutoff The cutoff level to filter insignificant p-values
#' 
#' @import stats
#' @examples 
#' data(test_data)
#' data(nonspecific_probes)
#' test_data_filtered <- filter_data(test_data)
#' test_data_ttest <- ttest_data(test_data_filtered, 1, 2, 3, 4, 1e-3)
#' 
#' @export

ttest_data <- function(x, cases_column_1, cases_column_n, controls_column_1, controls_column_n, ttest_cutoff) {
  x_ttest <- apply(x, 1, function(x) {t.test(x[cases_column_1:cases_column_n], x[controls_column_1:controls_column_n], "two.sided", var.equal = FALSE)$p.value})
  my_ttest_sorted <- sort(x_ttest, decreasing = FALSE)
  my_ttest_sorted_dtfm <- as.data.frame(my_ttest_sorted)
  my_ttest_candidate <- subset(my_ttest_sorted_dtfm, my_ttest_sorted_dtfm$my_ttest_sorted <= ttest_cutoff)
  return(my_ttest_candidate)
}