#' Applying delta beta analysis to calculate the difference between cases and controls
#' 
#' @description
#' This function calculates the delta beta value for the filtered probes. It calculates the difference in mean DNA methylation between cases and controls for each probe. Also, it selects probes with DNA methylation differences that are higher in cases than controls by a user specified meth_cutoff value and differences that are lower in cases than controls by the unmeth_cutoff value. In addition, the function provides the option to specify probes where the average beta value of the cases or controls is greater than a high_meth cutoff value or less than a low_meth cutoff value.
#' 
#' @param x the filtered 450k probes from filter_data function
#' @param cases_column_1 The first column (column number) for cases in the filtered dataset 
#' @param cases_column_n The last column (column number) for cases in the filtered dataset
#' @param controls_column_1 The first column (column number) for controls in the filtered dataset
#' @param controls_column_n The last column (column number) for controls in the filtered dataset
#' @param meth_cutoff The cutoff level for the methylation difference between cases and controls (cases minus controls)
#' @param unmeth_cutoff The cutoff level for the methylation difference between controls and cases (cases minus controls). Consequently, it requires a negative value.
#' @param high_meth The upper margin for the highly methylated probes
#' @param low_meth The lower margin for the low methylation
#' 
#' @examples 
#' data(test_data)
#' data(nonspecific_probes)
#' test_data_filtered <- filter_data(test_data)
#' test_data_delta_beta <- delta_beta_data(test_data_filtered, 1, 2, 3, 4, 0.5, -0.5, 0.94, 0.06)
#' 
#' @export


delta_beta_data <- function(x, cases_column_1, cases_column_n, controls_column_1, controls_column_n, meth_cutoff, unmeth_cutoff, high_meth, low_meth) {
  delta_beta_x <- rowMeans(as.matrix(x[,cases_column_1:cases_column_n])) - rowMeans(as.matrix(x[,controls_column_1:controls_column_n]))
  delta_beta_dtfm <- as.data.frame(delta_beta_x)
  delta_beta_meth_unmeth <- subset(delta_beta_dtfm, delta_beta_x >= meth_cutoff | delta_beta_x <= unmeth_cutoff)
  cases_mean <- rowMeans(as.matrix(x[,cases_column_1:cases_column_n]))
  cases_mean_dtfm <- as.data.frame(cases_mean)
  cases_meth_unmeth <- subset(cases_mean_dtfm, cases_mean >= high_meth | cases_mean <= low_meth)
  controls_mean <- rowMeans(as.matrix(x[,controls_column_1:controls_column_n]))
  controls_mean_dtfm <- as.data.frame(controls_mean)
  controls_meth_unmeth <- subset(controls_mean_dtfm, controls_mean >= high_meth | controls_mean <= low_meth)
  cases_controls <- merge(cases_meth_unmeth, controls_meth_unmeth, by = "row.names", all = TRUE)
  cases_control_rnames <- cases_controls[,1]
  cases_controls_m <- data.matrix(cases_controls[,2:ncol(cases_controls)])
  row.names(cases_controls_m) <- cases_control_rnames
  delta_cases_controls <- merge(delta_beta_meth_unmeth, cases_controls_m, by = "row.names")
  delta_cases_controls_rnames <- delta_cases_controls[,1]
  delta_cases_controls_m <- data.matrix(delta_cases_controls[,2])
  row.names(delta_cases_controls_m) <- delta_cases_controls_rnames
  return(delta_cases_controls_m)
}