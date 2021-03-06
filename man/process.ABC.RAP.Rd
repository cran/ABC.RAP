% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.ABC.RAP.R
\name{process.ABC.RAP}
\alias{process.ABC.RAP}
\title{An automated analysis applying all ABC.RAP functions in one script}
\usage{
process.ABC.RAP(x, cases_column_1, cases_column_n, controls_column_1,
  controls_column_n, ttest_cutoff, meth_cutoff, unmeth_cutoff, high_meth,
  low_meth)
}
\arguments{
\item{x}{The normalised beta values in a data matrix format, where conditions are arranged in columns and cg probes are arranged in rows.}

\item{cases_column_1}{The first column (column number) for cases in the filtered dataset}

\item{cases_column_n}{The last column (column number) for cases in the filtered dataset}

\item{controls_column_1}{The first column (column number) for controls in the filtered dataset}

\item{controls_column_n}{The last column (column number) for controls in the filtered dataset}

\item{ttest_cutoff}{The cutoff level to filter insignificant p-values}

\item{meth_cutoff}{The cutoff level for the methylation difference between cases and controls (cases minus controls)}

\item{unmeth_cutoff}{The cutoff level for the methylation difference between controls and cases (controls minus cases). Consequently, it requires a negative value.}

\item{high_meth}{The upper margin for the highly methylated probes}

\item{low_meth}{The lower margin for the low methylation}
}
\description{
This function processes the ABC.RAP workflow automatically
}
\examples{
\donttest{
data(test_data)
data(nonspecific_probes)
data(annotation_file)
process.ABC.RAP(test_data, 1, 2, 3, 4, 1e-3,  0.5, -0.5, 0.94, 0.06)
}
}

