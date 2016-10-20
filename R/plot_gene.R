#' Plotting and exporting methylation profile for candidate genes 
#' 
#' @description 
#' This function explores the DNA methylation profile for any gene. The function generates four plots: the top plots show the difference in DNA methylation between cases and controls (a bar chart of the delta beta values for all probes arranged from 5’ to 3’ positions and a plot showing the difference in mean DNA methylation between cases and controls). The bottom plots show the distribution of DNA methylation for each probe that interrogates a CpG site in the investigated gene, for cases (left) and controls (right), respectively. Also, an annotation table for the arranged probes is generated with the following columns: probe names, gene name, distance from TSS, mean methylation for cases, mean methylation for controls, delta beta values (cases minus controls), and t-test p.values.
#' 
#' @param x The filtered and annotated 450k probes 
#' @param b Gene name between quotation marks
#' @param cases_column_1 The first column (column number) for cases in the filtered dataset 
#' @param cases_column_n The last column (column number) for cases in the filtered dataset
#' @param controls_column_1 The first column (column number) for controls in the filtered dataset
#' @param controls_column_n The last column (column number) for controls in the filtered dataset
#' 
#' @examples  
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' test_data_filtered <- filter_data(test_data)
#' test_data_annotated <- annotate_data(test_data_filtered)
#' KLHL34 <- plot_gene(test_data_annotated, 'KLHL34', 1, 2, 3, 4)
#' 
#' @export


plot_gene <- function(x, b, cases_column_1, cases_column_n, controls_column_1, controls_column_n) {
  b_gene <- x[grep(b, x$Gene),]
  b_ordered <- b_gene[order(b_gene[,3]),]
  b_beta <- b_ordered[,6:ncol(b_ordered)]
  delta_beta <- rowMeans(as.matrix(b_beta[,cases_column_1:cases_column_n])) - rowMeans(as.matrix(b_beta[,controls_column_1:controls_column_n]))
  cases_betas <- b_beta[,cases_column_1:cases_column_n]
  controls_betas <- b_beta[,controls_column_1:controls_column_n]
  cases_mean <- rowMeans(as.matrix(b_beta[,cases_column_1:cases_column_n]))
  controls_mean <- rowMeans(as.matrix(b_beta[,controls_column_1:controls_column_n]))
  ttest_p.value <- apply(b_beta, 1, function(x) {t.test(x[cases_column_1:cases_column_n], x[controls_column_1:controls_column_n], "two.sided", var.equal = FALSE)$p.value})
  b_exported <- cbind(b_ordered[,c(1,2,4,5)], cases_mean, controls_mean, delta_beta, ttest_p.value)
  barplot_cols <- c("blue", "red")[(delta_beta > 0) + 1]
  par(mfrow=c(2,2))
  barplot(delta_beta, names.arg = b_beta$row.names, col = barplot_cols, main = paste("difference in DNA methylation of cases minus controls for", toString(b)), xlab = "450k probes arranged from 5' -> 3'", ylab = "delta beta values")
  plot(cases_mean, col = "red", ylim = c(0, 1),  main = paste("mean DNA methylation for",toString(b)), sub = "red circles = cases, blue triangles = controls", xlab = "450k probes arranged from 5' -> 3'", ylab = "beta values")
  points(controls_mean, pch = 24, col = "blue", bg = "blue")
  boxplot(t(cases_betas), col = "red", ylim = c(0, 1), main = paste(toString(b), ":DNA methylation for cases"), xlab = "450k probes arranged from 5' -> 3'", ylab = "beta values")
  boxplot(t(controls_betas), col = "blue", ylim = c(0,1), main = paste(toString(b), ":DNA methylation for controls"), xlab = "450k probes arranged from 5' -> 3'", ylab = "beta values")
  print(b_exported)
}