% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_igv.R
\name{runBasic}
\alias{runBasic}
\title{run igvR to plot SNP p-values over gene models}
\usage{
runBasic(
  chr = "20",
  start = 60000,
  end = 1100000,
  dn8like,
  use = "mlog10p",
  posvbl = "pos_b38",
  pvbl = "p"
)
}
\arguments{
\item{chr}{character(1)}

\item{start}{numeric(1)}

\item{end}{numeric(1)}

\item{dn8like}{data.frame}

\item{use}{character(1) defaults to "mlog10p"}

\item{posvbl}{character(1) what column in dn8like table holds genomic position}

\item{pvbl}{character(1) what column in dn8like table holds p-value}
}
\description{
run igvR to plot SNP p-values over gene models
}
\examples{
data("gtex_b38_lung_chr20_exc", package="tnt4dn8")
runBasic(dn8like=gtex_b38_lung_chr20_exc)
}
