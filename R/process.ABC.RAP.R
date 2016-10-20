#' An automated analysis applying all ABC.RAP functions in one script
#' 
#' @description 
#' This function processes the ABC.RAP workflow automatically
#' 
#' @param x The normalised beta values in a data matrix format, where conditions are arranged in columns and cg probes are arranged in rows.
#' @param cases_column_1 The first column (column number) for cases in the filtered dataset 
#' @param cases_column_n The last column (column number) for cases in the filtered dataset
#' @param controls_column_1 The first column (column number) for controls in the filtered dataset
#' @param controls_column_n The last column (column number) for controls in the filtered dataset
#' @param ttest_cutoff The cutoff level to filter insignificant p-values
#' @param meth_cutoff The cutoff level for the methylation difference between cases and controls (cases minus controls)
#' @param unmeth_cutoff The cutoff level for the methylation difference between controls and cases (controls minus cases). Consequently, it requires a negative value.
#' @param high_meth The upper margin for the highly methylated probes
#' @param low_meth The lower margin for the low methylation
#' 
#' @import grDevices
#' @examples 
#' \donttest{
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' process.ABC.RAP(test_data, 1, 2, 3, 4, 1e-3,  0.5, -0.5, 0.94, 0.06)
#' }
#' @export

process.ABC.RAP <- function(x, cases_column_1, cases_column_n, controls_column_1, controls_column_n, ttest_cutoff, meth_cutoff, unmeth_cutoff, high_meth, low_meth) {
  x_filtered <- filter_data(x)
  x_annotated <- annotate_data(x_filtered)
  x_ttest <- ttest_data(x_filtered, cases_column_1, cases_column_n, controls_column_1, controls_column_n, ttest_cutoff)
  x_delta <- delta_beta_data(x_filtered, cases_column_1, cases_column_n, controls_column_1, controls_column_n, meth_cutoff, unmeth_cutoff, high_meth, low_meth)
  x_overlap <- overlap_data(x_ttest, x_delta)
  x_CpG_hits <- CpG_hits(x_overlap)
  pdf(file = "process.ABC.RAP.plots.pdf", width = 11)
  sink("process.ABC.RAP.tables.txt")
  for(i in x_CpG_hits$Var1) {
    plot_gene(x_annotated, i, cases_column_1, cases_column_n, controls_column_1, controls_column_n)
  }
  dev.off()
  sink()
}