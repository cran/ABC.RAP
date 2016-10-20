utils::globalVariables('barplot')
#' Plotting highly different and significant probes annotated by their corresponding gene names
#' 
#' @description 
#' This function plots the potential candidate genes for which multiple CpG sites show significant difference.
#' @param x Results from the overlap_data function
#' @import graphics
#' @examples 
#' \donttest{
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' test_data_filtered <- filter_data(test_data)
#' test_data_ttest <- ttest_data(test_data_filtered, 1, 2, 3, 4, 1e-3)
#' test_data_delta_beta <- delta_beta_data(test_data_filtered, 1, 2, 3, 4, 0.5, -0.5, 0.94, 0.06)
#' test_overlapped_data <- overlap_data(test_data_ttest, test_data_delta_beta)
#' plot_candidate_genes(test_overlapped_data)
#' }
#' @export


plot_candidate_genes <- function(x) {
  ttest_delta_table <- data.frame(table(x$Gene))
  ttest_delta_multiple_hits <- x[x$Gene %in% ttest_delta_table$Var1[ttest_delta_table$Freq>1],]
  ttest_delta_multiple_hits_ordered <- ttest_delta_multiple_hits[order(ttest_delta_multiple_hits[,2]),]
  par(mfrow = c(1,1))
  cols <- c("blue", "red")[(ttest_delta_multiple_hits_ordered$V1 > 0) + 1]
  barplot(ttest_delta_multiple_hits_ordered$V1, names.arg = ttest_delta_multiple_hits_ordered$Gene, col = cols, main = "Genes with multiple CpG sites", ylab = "delta beta DNA methylation values")
}