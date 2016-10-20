#' Filtering DNA methylation 450k non_specific probes
#' 
#' @description 
#' This function filters the reported nonspecific probes, and also filters probes that interrogate SNPs of minor allele frequency (MAF) > 0.1. A list of nonspecific probes was obtained from Chen et al (2013) supplementary files.
#' @param x The normalised beta values in a data matrix format, where conditions are arranged in columns and cg probes are arranged in rows. 
#'
#' @references Chen YA, Lemire M, Choufani S, et al. Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. Epigenetics 2013;8:203-9. 
#' @examples
#' data(test_data)
#' data(nonspecific_probes)
#' test_data_filtered <- filter_data(test_data)
#'
#' @export


filter_data <- function(x) {
  x_filtered <- x[!(x[,1] %in% nonspecific_probes[,1]),]
  x_rnames <- x_filtered[,1]
  x_m <- data.matrix(x_filtered[,2:ncol(x_filtered)])
  row.names(x_m) <- x_rnames
  return(x_m)
}