#' Overlapping Student's t-test and delta beta results
#' 
#' @description 
#' This function overlaps the results from both Student’s t-test and delta beta analyses to identify probes (CpG sites) that are highly and significantly different between cases and controls.
#'
#' @param x Results from t-test or delta beta analyses 
#' @param y Results from t-test or delta beta analyses 
#' @examples 
#' \donttest{
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' test_data_filtered <- filter_data(test_data)
#' test_data_ttest <- ttest_data(test_data_filtered, 1, 2, 3, 4, 1e-3)
#' test_data_delta_beta <- delta_beta_data(test_data_filtered, 1, 2, 3, 4, 0.5, -0.5, 0.94, 0.06)
#' test_overlapped_data <- overlap_data(test_data_ttest, test_data_delta_beta)
#' }
#' @export


overlap_data <- function(x, y) {
  ttest_delta_candidate <- merge(x, y, by = "row.names")
  ttest_delta_rnames <- ttest_delta_candidate[,1]
  ttest_delta_m <- data.matrix(ttest_delta_candidate[,2:ncol(ttest_delta_candidate)])
  row.names(ttest_delta_m) <- ttest_delta_rnames
  annotation_file_rnames <- annotation_file[,1]
  annotation_f <- data.frame(annotation_file[,2:ncol(annotation_file)])
  row.names(annotation_f) <- annotation_file_rnames
  ttest_delta_annotated <- merge.data.frame(annotation_f, ttest_delta_m, by = "row.names")
  ttest_delta_arranged <- ttest_delta_annotated[order(ttest_delta_annotated[,4]),]
  return(ttest_delta_arranged)
}