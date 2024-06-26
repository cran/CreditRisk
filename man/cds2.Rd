% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cds2.R
\name{cds2}
\alias{cds2}
\title{Calculate Credit Default Swap rates}
\usage{
cds2(t, Tj, tr, r, tint, int, R = 0.005, ...)
}
\arguments{
\item{t}{premium timetable.}

\item{Tj}{CDS maturities.}

\item{tr}{interest rates timetable.}

\item{r}{spot interest rates.}

\item{tint}{intensity timetable.}

\item{int}{default intensities timetable.}

\item{R}{constant premium payment.}

\item{...}{further arguments on \code{cds}.}
}
\value{
An object of class \code{data.frame} that contains the quantities calculated by \code{cds}
on Tj timetable.
}
\description{
Calculate CDS rates starting from default intensities
}
\details{
The function \code{cds2} is based on \code{cds} but allows a more fine controll on maturities
and on discretization of \code{r} and \code{int}. In particular input \code{(t, tr, tint)}
can be of different length thanks to the function \link[stats]{approx}.
}
\examples{
cds2(t = c(1:20),Tj = c(1:20), tr = c(1:20), r = seq(0.01,0.06, len =20),
tint = c(1:20), int= seq(0.01,0.06, len =20))

}
\references{
David Lando (2004) Credit Risk Modeling.

Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
Counterparty Credit Risk, Collateral and Funding.
With Pricing Cases for All Asset Classes
}
