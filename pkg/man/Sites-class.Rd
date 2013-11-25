\name{Sites-class}
\Rdversion{1.1}
\docType{class}

\alias{sites}

\alias{Sites-class}

\alias{$,Sites-methods}

\alias{rbind,Sites-methods}

\alias{sites,Sites-methods}

\alias{sites-method}
\alias{sites,character-method}
\alias{sites,data.frame-method}
\alias{sites,matrix-method}
\alias{sites,Sites-method}
\alias{sites,rbind-method}

\title{Class \code{"Sites"}}

\description{
  Data representation of sites data in \emph{long format} where
  each observation on a particular plot (environmental variables,
  stand parameters) is represented by one row.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Sites", ...)},
  where \code{...} expands to slot \code{'data'}. Using the
  \code{sites} method with appropriate objects is also fast and convenient
  (see \sQuote{Methods}).

  The \code{data} slot holds a \code{"data.frame"} with three columns,
  all of type \code{character} and with names \code{plot}, \code{variable}
  and \code{value}. The meaning of these columns is as detailed below:

  \itemize{
    \item \code{plot} contains the plot identifier associated with each
      observation.
    \item \code{variable} contains a name for the measured variable.
    \item \code{value} holds the value of a particular variable.
  }
}

\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"}
      with three columns named
      \code{'plot'}, \code{'variable'} and \code{'value'}. }
  }
}

\section{Methods}{
  \describe{
    \item{ sites }{
      \code{signature(obj = "Sites")}:
      Accessor function (slot \code{obj@data}).

      \code{signature(obj = "data.frame")}:
      Promote \code{data.frame} to object of class \code{"Sites"}

      \code{signature(obj = "matrix")}:
      Promote \code{matrix} to object of class \code{"Sites"}
      Note, \code{obj} does not need to have \code{dimnames}.

      \code{signature(obj = "character")}:
        Supply a file to be read by \code{read.csv} and promote that object to
        class \code{"Sites"}. Remember to pass appropriate arguments to
        \code{\link{read.csv}} that meet the specifications of your file.
    }
    \item{ rbind }{
      \code{signature(x = "Sites")}:
        \code{rbind} method to combine two or more objects.
    }
  }
}

\author{ Roland Kaiser }

\seealso{
  \code{\link{stack.sites}}
  \code{\linkS4class{Species}}
  \code{\linkS4class{Vegsoup}}
}

\examples{
require(vegsoup)
showClass("Sites")

showClass("Species")

data(testdata)

new("Sites", data = sts)

class(sts)
sites(sts)
obj <- sites(as.matrix(sts))

# get slot data from object
class(sites(sts))

# combine objects
data(bigtestdata)

obj2 <- sites(sts.big)

obj <- rbind(obj, obj2)

}

\keyword{ classes }