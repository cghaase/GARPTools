% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commissionGARP.R
\name{commissionGARP}
\alias{commissionGARP}
\title{Calculates commission metrics of model outputs}
\usage{
commissionGARP(n, x, points, file.path)
}
\arguments{
\item{n}{a numeric value specifying the number of models in the best subset outputted by DesktopGARP}

\item{x}{a raster object of the summated raster of the best models output by DesktopGARP}

\item{points}{a spatial object of presence data to use for testing locations}

\item{file.path}{file directory of folder containing best model rasters output by DeskstopGARP (usually within athe "GARPruns folder")}
}
\value{
A list object containing:
 \itemize{
  \item \code{Average.Commission} the average commission of the model subset
  \item \code{Total.Commission} the total commission of the model subset
 }
}
\description{
\code{commissionGARP} Calculates the commission, or the proportion of the landscape predicted as present, of the best models output by DesktopGARP.
}
\details{
The discrete cutpoints (\code{n}) represent the number of models that agree on a predicted presence location.

The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.

The shapefile \code{points} should presence locations that were not used by DesktopGARP for model training and those output by \code{\link{splitData}}.
}
\examples{
\dontrun{
commissionGARP(n = 10, x = sum.raster, points = test_1, file.path = "C/GARPruns/BestSubsets/")
}

}
\seealso{
\code{\link{omissionGARP}}
}
