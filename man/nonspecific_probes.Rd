\name{nonspecific_probes}
\alias{nonspecific_probes}
\docType{data}
\title{450k DNA methylation non specific probes}
\description{data frame of the non specific probes that need to be filtered out from 450k array datasets}
\usage{data("nonspecific_probes")}
\format{A data frame}
\details{These non specific probes interrogates SNPs with mean allelic frequency (MAF) > 0.1, and also those that don't align uniquely on the genome. The list of nonspecific probes was obtained from Chen et al (2013) supplementary files}
\keyword{datasets}
\references{Chen YA, Lemire M, Choufani S, et al. Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. Epigenetics 2013;8;203-9}
