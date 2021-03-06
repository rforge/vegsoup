\encoding{UTF-8}

\name{Fidelity}
\docType{methods}

\alias{fidelity}
\alias{fidelity,VegsoupPartition-method}

\alias{FisherTest}
\alias{FisherTest,VegsoupPartition-method}
\alias{FisherTest,VegsoupPartitionFidelity-method}

\alias{Phi}
\alias{Phi,VegsoupPartition-method}

\alias{Indval}
\alias{Indval,VegsoupPartition-method}

\alias{getStat}
\alias{getStat,VegsoupPartitionFidelity-method}

\title{ Calculate Fidelity Values }

\description{
  Method \code{fidelity} calculates fidelity values for objects inheriting
  \code{'VegsoupPartition'} and returns an extended object of class
  \code{\linkS4class{VegsoupPartitionFidelity}}. Function \code{fidelity} is
  based on the implementation of function \code{\link[indicspecies]{strassoc}}
  in \pkg{indicspecies} and offers all it's methods plus some more
  (see \sQuote{Details}).

  Functions \code{Phi} and \code{Indval} calculates the \emph{phi} and \emph{Indval}
  statistics according to \code{\link[optpart:phi]{phi}} in package
  \pkg{optpart}, respectively, \code{\link[labdsv:indval]{indval}} of package
  \pkg{labdsv}. Function \code{FisherTest} is a slightly faster alternative to
  \code{fidelity(obj, method = "fisher")}. \code{Phi} and \code{Indval} are
  also faster implementations compared to \code{fidelity} without option
  \code{fast = TRUE}. All three methods do not return an \code{'VegsoupPartitionFidelity'}
  object, but the requested statistic as a \code{"data.frame"}! Take \code{sqrt(Indval(x))}
  to obtain comparable results with \code{fidelity(x, method = "IndVal.g")}.
}

\usage{
\S4method{fidelity}{VegsoupPartition}(obj, method = "r.g", group = NULL, nboot = 0, alpha = 0.05, c = 1,
		 alternative = "greater", fast = FALSE, verbose = TRUE, ...)

\S4method{FisherTest}{VegsoupPartition}(obj, alternative = "greater")

\S4method{FisherTest}{VegsoupPartitionFidelity}(obj)

\S4method{Phi}{VegsoupPartition}(obj)

\S4method{Indval}{VegsoupPartition}(obj, ...)

\S4method{getStat}{VegsoupPartitionFidelity}(obj)
}

\arguments{
  \item{ obj }{
	\code{'VegsoupPartition'} object.
   }
  \item{ method }{
	indicator measure.
  }
  \item{ group }{
	grouping variable.
  }
  \item{ nboot }{
	number of boot straps, defaults to.
  }
  \item{ alpha }{
	significance.
  }
  \item{ c }{
	Forgot what that this is?
  }
  \item{ fast }{
	accelerate computations using package \code{parallel}.
  }
  \item{ verbose }{
	display computation status bar.
  }
  \item{ alternative }{
	alternative hypothesis. For \code{fidelity} takes effect only for \code{method = "Fisher"}
  }
  \item{ \dots }{
	For \code{fidelity} additional arguments for \code{method = "fisher.s"}
	and for \code{Indval} arguments passed to \code{\link[labdsv:indval]{indval}}.
  }
}

\details{
Methods currently implemented for method \code{fidelity} include the following:

  \itemize{
	 \item \code{r:} phi (point-biserial correlation coefficient)
	 \item \code{r.g:} phi, group equalised
	 \item \code{cos:} cosine (Ochiai index)
	 \item \code{cos.g:} cosine (Ochiai index), group equalised
	 \item \code{r.ind:} phi individual based
	 \item \code{r.ind.g:} phi individual based, group equalised
	 \item \code{s.ind:} square root of Indval A times Indval B
	 \item \code{s.ind.g:} square root of Indval A.g * Indval B, group equalised
	 \item \code{IndVal.g:} Indval, group equalised
	 \item \code{IndVal:} Indval
	 \item \code{A.g:} Indval A, group equalised
	 \item \code{A:} Indval A
	 \item \code{B:} Indval B
	 \item \code{u.hyp:} Bruelheide's corrected u value (Bruelheide 2000)
	 \item \code{u.binB:} Bruelheide's corrected u value uBinB (Bruelheide 2000)
	 \item \code{u.binA:} Bruelheide's corrected u value uBinA (Bruelheide 2000)
	 \item \code{g:} G statistic
	 \item \code{chi:} chi-square statistic with Yates correction
	 \item \code{Fisher:} Fisher test (Chytry et al. 2002)
	 \item \code{CR:} Willner's constancy ratio (Willner 2006)
	 \item \code{TCR:} Willner's total cover ratio (Willner 2006)
  }
}
\value{
For \code{fidelity} an object of class \code{'VegsoupPartitionFidelity'}.

For \code{FisherTest} and \code{Phi} a \code{'matrix'} of values with species as rows and partitions as columns. Note, for class \code{\linkS4class{VegsoupPartitionFidelity}} \code{FisherTest} just extracts what was already calculated.
}
\references{
Willner, W., Tichy, L., and Chytry, M. (2009). Effects of different fidelity measures and contexts on the determination of diagnostic species. Journal of Vegetation Science, 20(1):130-137.

Willner, W. (2006). The association concept revisited. Phytocoenologia, 36(1):67-76.

Tichy, L. and M. Chytry. 2006. Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science 17:809-818.

Chytry, M., Tichy, L., Holt, J., and Botta-Dukat, Z. (2002). Determination of diagnostic species with statistical fidelity measures. Journal of Vegetation Science, 13(1):79-90.

De \enc{Cáceres}{Caceres}, M., Font, X., and Oliva, F. (2008). Assessing species diagnostic value in large data sets: A comparison between phi-coefficient and ochiai index. Journal of Vegetation Science, 19(6):779-788.

De \enc{Cáceres}{Caceres}, M. and Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12):3566-3574.

De \enc{Cáceres}{Caceres}, M., Legendre, P., and Moretti, M. (2010). Improving indicator species analysis by combining groups of sites. Oikos, 119(10):1674-1684.

Bruelheide, H. (2000). A new measure of fidelity and its application to defining species groups. Journal of Vegetation Science, 11(2):167-178.
}
\author{
Roland Kaiser, based on code from Miquel De \enc{Cáceres}{Caceres} Ainsa taken from function \code{\link[indicspecies:strassoc]{strassoc}} in package \pkg{indicspecies}.
}
\seealso{
\code{\linkS4class{VegsoupPartitionFidelity}}
}
\note{
  Calculations can be computationally intensive -- be patient.
}
\examples{
require(vegsoup)

data(windsfeld)

x <- VegsoupPartition(windsfeld, k = 2)
f <- fidelity(x)

head(getStat(f))

hist(f)
}
\keyword{ methods }
\keyword{ multivariate }

