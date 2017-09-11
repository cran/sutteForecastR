\name{alpha.sutte}
\alias{alpha.sutte}
\title{Forecasting Data using Alpha-Sutte Indicator}
\description{
The alpha-Sutte indicator (alpha-Sutte) was originally from developed of Sutte indicator. Sutte indicator can using to predict the movement of stocks. As the development of science, then Sutte indicator developed to predict not only the movement of stocks but also can forecast data on financial, insurance, and others time series data.
}
\usage{
alpha.sutte(x)
}
\arguments{
\item{x}{ a time series data}
}
\author{
Ansari Saleh Ahmar <ansarisaleh@unm.ac.id>
}
\references{
Ahmar, A.S. (2017) <doi:10.17605/osf.io/rknsv>

Ahmar, A. S., Rahman, A., & Mulbar, U. (2017) <doi:10.17605/osf.io/n68yq>
}
\examples{
x <- c(94.77, 96.23, 98.12, 99.09, 100.04, 100.12, 99.93, 100.09,
101.44, 102.38, 103.68,104.12, 104.81, 105.35, 106.36, 106.89,107.35,
107.21, 107.72, 108.54, 109.57, 112.03)

alpha.sutte(x)
}
