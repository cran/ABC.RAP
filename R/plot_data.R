#' Overview description of the DNA methylation pattern for cases and controls
#' 
#' @description 
#' This function produces four distribution plots that summarise the DNA methylation patterns for cases (top left) and controls (top right). The top two histograms show the pattern of mean DNA methylation levels for cases and controls. The bottom two plots show the difference in DNA methylation between cases and controls (a boxplot comparing methylation profile for cases and controls, and a delta beta plot describing the methylation difference between cases and controls). The function also provides summary statistics for the delta beta analysis that can be used to select cutoff values for the delta_beta_data function.
#' 
#' @param x The filtered 450k probes from filter_data() function
#' @param cases_column_1 The first column (column number) for cases in the filtered dataset 
#' @param cases_column_n The last column (column number) for cases in the filtered dataset
#' @param controls_column_1 The first column (column number) for controls in the filtered dataset
#' @param controls_column_n The last column (column number) for controls in the filtered dataset
#' 
#' @examples 
#' data(test_data)
#' data(nonspecific_probes)
#' test_data_filtered <- filter_data(test_data)
#' plot_data(test_data_filtered, 1, 2, 3, 4)
#' 
#' @export


plot_data <- function(x, cases_column_1, cases_column_n, controls_column_1, controls_column_n) {
  cases <- x[,cases_column_1:cases_column_n]
  controls <- x[,controls_column_1:controls_column_n]
  cases_mean <- rowMeans(as.matrix(x[,cases_column_1:cases_column_n]))
  controls_mean <- rowMeans(as.matrix(x[,controls_column_1:controls_column_n]))
  delta_beta <- rowMeans(as.matrix(x[,cases_column_1:cases_column_n])) - rowMeans(as.matrix(x[,controls_column_1:controls_column_n]))
  means <- cbind(cases_mean, controls_mean)
  Hist_cal <- hist(delta_beta, plot = FALSE)
  hist_categ <- cut(Hist_cal$breaks, c(-Inf, -0.00001, 0.0, Inf))
  par(mfrow=c(2,2))
  hist(cases, col = "red", xlab = "DNA methylation")
  hist(controls, col = "blue", xlab = "DNA methylation")
  boxplot(means, col = c("red","blue"), ylim = c(0,1), main = "DNA methylation of cases and controls")
  plot(Hist_cal, col = c("blue", "red", "red")[hist_categ], main = "Difference of DNA methylation between cases and controls", xlab = "delta beta values of cases minus controls")
  print(summary(delta_beta))
}