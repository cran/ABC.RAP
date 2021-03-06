% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_data.R
\name{plot_data}
\alias{plot_data}
\title{Overview description of the DNA methylation pattern for cases and controls}
\usage{
plot_data(x, cases_column_1, cases_column_n, controls_column_1,
  controls_column_n)
}
\arguments{
\item{x}{The filtered 450k probes from filter_data() function}

\item{cases_column_1}{The first column (column number) for cases in the filtered dataset}

\item{cases_column_n}{The last column (column number) for cases in the filtered dataset}

\item{controls_column_1}{The first column (column number) for controls in the filtered dataset}

\item{controls_column_n}{The last column (column number) for controls in the filtered dataset}
}
\description{
This function produces four distribution plots that summarise the DNA methylation patterns for cases (top left) and controls (top right). The top two histograms show the pattern of mean DNA methylation levels for cases and controls. The bottom two plots show the difference in DNA methylation between cases and controls (a boxplot comparing methylation profile for cases and controls, and a delta beta plot describing the methylation difference between cases and controls). The function also provides summary statistics for the delta beta analysis that can be used to select cutoff values for the delta_beta_data function.
}
\examples{
data(test_data)
data(nonspecific_probes)
test_data_filtered <- filter_data(test_data)
plot_data(test_data_filtered, 1, 2, 3, 4)

}

