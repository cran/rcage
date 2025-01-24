% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wendland.R
\name{wendland}
\alias{wendland}
\title{Wendland Basis Functions}
\usage{
wendland(crd, knots, w = NULL, ..., longlat = TRUE)
}
\arguments{
\item{crd}{A matrix object. The (x,y) coordinates of the reference
points \{nCrd x 2\}.}

\item{knots}{A matrix object. The (x,y) coordinates of the knots \{r x 2\}.}

\item{w}{A numeric object. The positive scaling factor (bandwidth).}

\item{...}{ignored. Included only to require naming of inputs that follow.}

\item{longlat}{A logical object. If FALSE, Euclidean
distance is calculated; if TRUE, Great Circle distance is calculated.
See ?sp::spDists for more information.}
}
\value{
A matrix of Wendland functions evaluated at all combinations of
 crd and knots \{ nCrd x r \}.
}
\description{
Provides an implementation of the Wendland basis functions defined as
\deqn{ \Psi_{j}(s) =
 \{ 1 - d_{j}(s)\}^6 \{35 d_{j}(s)^2 + 18 d_j(s) + 3\}/3 \mathrm{I}( 0 \leq d_{j} \leq 1  ),}
where 
\deqn{ d_{j}(s) = ||s - c_j||/w.}
}
\details{
Distances between reference coordinates and knots are obtained using 
  sp::spDists().
}
\examples{

  data(countyExample)

  wendland(crd = sp::coordinates(county), knots = knots)

}
\references{
Wendland, H. (1998). Error estimates for interpolation by
  compactly supported radial basis functions of minimal degree. Journal
  of Approximation Theory, 93,258-272. <doi:10.1006/jath.1997.3137>.
}
