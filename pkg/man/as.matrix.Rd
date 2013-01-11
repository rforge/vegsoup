\name{as.matrix}
\alias{as.logical}
\alias{as.numeric}
\alias{as.character}
\alias{as.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve Species Matrix with Specified Storage Mode
}
\description{
The functions retrieve a species matrix of mode according the methods name and \code{as.matrix} returns corresponding types of matrices but specified by its argument \code{mode}.
}
\usage{
as.logical(x, mode)

as.numeric(x, mode)

as.character(x, mode)

as.matrix(x, typeof, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object inheriting from class \code{VegsoupData}.
}
  \item{\dots}{
If \code{mode = "R"} return the transposed matrix. Defaults to \code{"Q"}, no transposition. See \sQuote{Details}. 
}
  \item{typeof}{
A character string giving the desired storage mode or type of the object. Can be one of \code{"character"}, \code{"numeric"} or \code{"logical"}.
}
  \item{\dots}{
Arguments passed to \code{as.logical} etc. to supply a value for the argument \code{mode}.
}
}
\details{
For any object inheriting from class \code{VegsoupData} the methods return an object of class \code{"matrix"} with dimnames giving plot and species names. The latter are strings where any occurrence in a layer is labeled using the layer name and the species abbreviations, collapased with the \code{"@"} gylph (compound abbreviated taxa names, see \code{\link{DecomposeNames}} and \code{\link{Taxonomy}}). The returned object is suitable as input to functions in other packages that require an object of class \code{matrix}. Allmost all of these functions internally coerce objects of class \code{"data.frame"} to class \code{"matrix"}.

\code{as.numeric} returns a matrix of mode \code{"numeric"} with species abundances coded by \code{AbundanceScale(obj)\$lims}. Absences are zeros. Equals \code{as.matrix(x, typeof = "numeric")}.

\code{as.character} returns a matrix of mode \code{"character"} with species abundances in the data set's original scale \code{AbundanceScale(obj)\$codes}. Absences are coded as zeros! Equals \code{as.matrix(x, typeof = "character")}.

\code{as.logical} returns a matrix of mode \code{"numeric"} with zeros coding for species absences and ones for presences. Equals \code{as.matrix(obj, mode = "logical")}.

Note, \code{as.logical} and \code{as.matrix(x, "logical")} do not return a boolean or logical matrix composed of \code{TRUE} and \code{FALSE} but use ones to code presences (\code{TRUE}) and zeros (\code{FALSE}) for absences. The \code{mode} (\code{storage.mode}) of the returned matrix is of type \code{"numeric"} (\code{"integer"}).

As defined above, the some how misleading named argument \code{mode} refers to the type of analysis intended on the returned object (recall \code{"typeof"} determines \code{mode} (\code{storage.mode}) of the returned matrix). The argument \code{mode} can be passed to all methods and has two possible values. If \code{mode = "Q"} pairs of objects are compared and the analysis is said to be in \emph{Q mode}. This is what users usally expect. When pairs of descriptors are compared, the analysis is said to be in \emph{Q mode} (Legendre & Legendre 2012). For example, \code{as.matrix(x, "logical", "R")} returns a transposed presence/absence matrix. Of cource, this can also be written as \code{t{as.matrix(x, "logical")}}.

}
\value{
A matrix of mode \code{"numeric"} or \code{"character"} (see \sQuote{Details}).
}
\note{
If the object has non \code{NULL} values in its data standardization slot the method defined by \code{decostand(obj)} will automatically be applied!
}
\author{
Roland Kaiser
}
\references{
Legendre, P. and Legendre, L. (2012). Numerical ecology. Developments in environmental modelling, 24.
}
\seealso{
\code{\link{Vegsoup}}, \code{\link{VegsoupData}}, \code{\link{decostand}} and \code{\link{AbundanceScale}}
}
\examples{
library(vegsoup)
data(testdata)
dta <- VegsoupData(Vegsoup(species, sites, taxonomy,
	scale = list(scale = "Braun-Blanquet")))

#	presence/absence representation
dta.pa <- as.logical(dta)
#	the same results can be achieved using the classes decostand property
decostand(dta) <- "pa"

#	note, that the results of as.numeric() and as.logical()
#	are the same if presence/absence (transformation) is applied
dta.pa2 <- as.numeric(dta)

#	for a true boolean matrix change mode
mode(dta.pa) <- "logical"
head(dta.pa)
all.equal(dta.pa, dta.pa2, check.attributes = FALSE)

#	using as.matrix and argument mode
as.matrix(dta, "character", "q")


require(vegan)

pro <- procrustes(as.logical(dta), as.numeric(dta))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
