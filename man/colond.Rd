% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colond.R
\docType{data}
\name{colond}
\alias{colond}
\title{Dataset of colon cancer deaths.}
\format{A data frame with 888 rows and 13 variables:
\describe{
  \item{time}{Survival time (in days).}
  \item{status}{Survival time status: Observed (\code{TRUE}), Right-censored (\code{FALSE}).}
  \item{rx}{Treatment type: Control (\code{obs}), Levamisole (\code{Lev}), Levamisole + Fluorouracil (\code{Lev+5FU}).}
  \item{sex}{Sex of patient (\code{female}/\code{male}).}
  \item{age}{Age of patient (in years).}
  \item{obstruct}{Obstruction of colon by tumor (T/F).}
  \item{perfor}{Perforation of colon (T/F).}
  \item{adhere}{Adherence of tumor to nearby organs (T/F).}
  \item{nodes}{Number of lymph nodes with detectable cancer (integer betwwen 0-33).}
  \item{differ}{Differentiation index of tumor: \code{well}, moderate (\code{moder}), \code{poor}.}
  \item{extent}{Extent of local spread: submucosa (\code{submuc}), muscle, serosa, contiguous structures (\code{cstruct}).}
  \item{surg}{Time from surgery to registration (\code{short}, \code{long}).}
  \item{node4}{More than 4 positive lymph nodes (T/F).}
}}
\source{
Subset of the \code{colon} dataset from the \pkg{survival} package.
}
\usage{
colond
}
\description{
Survival times of patients in the chemotherapy colon cancer study of Laurie (1989).  This is a subset of the \code{colon} dataset from the \pkg{survival} package.
}
\references{
Laurie, J.A., Moertel, C.G., Fleming, T.R., Wieand, H.S., Leigh, J.E., Rubin, J., McCormack, G.W., Gerstner, J.B., Krook, J.E., and Malliard, J. "Surgical adjuvant therapy of large-bowel carcinoma: An evaluation of levamisole and the combination of levamisole and fluorouracil: The North Central Cancer Treatment Group and the Mayo Clinic".  \emph{Journal of Clinical Oncology} 7:10 (1989): 1447-1456.  \url{https://www.ncbi.nlm.nih.gov/pubmed/2778478}
}
\keyword{datasets}
