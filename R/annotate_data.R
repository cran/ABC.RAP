utils::globalVariables(c('annotation_file', 'test_data', 'nonspecific_probes'))
#' Annotating the filtered probes
#' 
#' @description
#' This function annotates each filtered probe with gene name, chromosome number, probe location, distance from transcription start site (TSS), and relation to CpG islands. The annotation file is based on "UCSC platform" annotation format and was obtained from Illumina GPL13534_HumanMethylation450_15017482_v1.1 file (BS0010894-AQP_content.bpm). 
#' @import utils
#' @param x the filtered probes from filter_data 
#' @examples 
#' data(test_data)
#' data(nonspecific_probes)
#' data(annotation_file)
#' test_data_filtered <- filter_data(test_data)
#' test_data_annotated <- annotate_data(test_data_filtered)
#' 
#' @export

annotate_data <- function(x) {
  annotation_file_rnames <- annotation_file[,1]
  annotation_f <- data.frame(annotation_file[,2:ncol(annotation_file)])
  row.names(annotation_f) <- annotation_file_rnames
  x_merged <- merge.data.frame(annotation_f, x, by = "row.names")
  x_rnmaes <- x_merged[,1]
  x_annotated <- data.frame(x_merged[,2:ncol(x_merged)])
  row.names(x_annotated) <- x_rnmaes
  return(x_annotated)
}