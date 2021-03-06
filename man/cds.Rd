% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cds.R
\name{cds}
\alias{cds}
\title{Calculates Credit Default Swap rates}
\usage{
cds(t, int, r, R = 0.005, RR = 0.4, simplified = FALSE)
}
\arguments{
\item{t}{premium timetable.}

\item{int}{deterministic default intensities vector.}

\item{r}{spot interest rates.}

\item{R}{constant premium payments, value that the buyer pays in each \eqn{t_i}.}

\item{RR}{recovery rate on the underline bond, default value is 40\%.}

\item{simplified}{logic argument. If FALSE calculates the CDS rates using the
semplified version of calculations, if TRUE use the complete version.}
}
\value{
\code{cds} returns an object of class \code{data.frame} with columns, for esch date
\eqn{t_i} the value of survival probability, the premium and protection leg, CDS rate
and CDS price.
}
\description{
Calculates CDS rates starting form default intensities.
}
\details{
\itemize{
\item  Premium timetable is \eqn{t_i;  i=1,...,T}. The vector starts from
\eqn{t_1\le 1}, i.e. the first premium is payed at a year fraction in the possibility that
 the bond is not yet defaulted. Since premium are a postponed payment (unlike usual insurance
  contracts).
\item Intensities timetable have domains \eqn{\gamma_i;  i=t_1,...,T}.
\item  spot interest rates of bond have domain \eqn{r_i; i=t_1,...,T}. The function transforms
 spot rates in forward rates. If we specify that we want to calculate CDS rates with the
 simplified alghoritm, in each period, the amount of the constant premium payment
 is expressed by: \deqn{\pi^{pb}=\sum_{i=1}^Tp(0,i)S(0,i)\alpha_i}
 and the amount of protection, assuming a recovery rate \eqn{\delta}, is:
 \deqn{\pi^{ps}=(1-\delta)\sum_{i=1}^Tp(0,i)\hat{Q}(\tau=i)\alpha_i}
 If we want to calculate same quantities with the complete version, that evaluate premium in the continous,
 the value of the premium leg is calculated as:
 \deqn{\pi^{pb}(0,1)=-\int_{T_a}^{T_b}P(0,t)\cdot(t-T_{\beta(t)-1}) d_t Q
 (\tau\geq t)+\sum_{i=a+1}^bP(0,T_i)\cdot\alpha_i * Q(\tau\geq T_i)} and the protection leg
 as:
 \deqn{\pi_{a,b}^{ps}(1):=-\int_{t=T_a}^{T_b}P(0,t)d*Q(\tau\geq t)}
 In both versions the forward rates and intensities are supposed as costant stepwise functions
  with discontinuity in \eqn{t_i}
}
}
\examples{
cds(t = seq(0.5, 10, by = 0.5), int = seq(.01, 0.05, len = 20),
r = seq(0,0.02, len=20), R = 0.005, RR = 0.4, simplified = FALSE)

}
\references{
David Lando  (2004) Credit risk modeling.

Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
Counterparty Credit Risk, Collateral and Funding.
With Pricing Cases for All Asset Classes
}
