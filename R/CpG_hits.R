#' Identifying genes for which multiple CpG sites show significant methylation difference
#' 
#' @description 
#' This function calculates the number of significantly different CpG sites between cases and controls for each gene and produces a frequency table with genes that have more than one CpG site. 
#' 
#' @param x Results from the overlap_data function
#' @examples
#' \donttest{
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' test_data_filtered <- filter_data(test_data)
#' test_data_ttest <- ttest_data(test_data_filtered, 1, 2, 3, 4, 1e-3)
#' test_data_delta_beta <- delta_beta_data(test_data_filtered, 1, 2, 3, 4, 0.5, -0.5, 0.94, 0.06)
#' test_overlapped_data <- overlap_data(test_data_ttest, test_data_delta_beta)
#' test_CpG_hits <- CpG_hits(test_overlapped_data)
#' }
#' @export

CpG_hits <- function(x) {
  ttest_delta_table <- data.frame(table(x$Gene))
  ttest_delta_multiple_hits <- x[x$Gene %in% ttest_delta_table$Var1[ttest_delta_table$Freq>1],]
  multiple_hits <- ttest_delta_table[ttest_delta_table$Freq > 1,]
  return(multiple_hits)
}

