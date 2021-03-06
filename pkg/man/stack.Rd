\encoding{UTF-8}

\name{stack}

\alias{stack}
\alias{stackSpecies}
\alias{stackSites}
\alias{stackCoordinates}

\title{ Stack Species Matrix, Sites Data or Spatial Coordinates }

\description{
  All functions reshape a data frame between \emph{wide format} with
  measurements in separate columns to \emph{long format} with the measurements
  in separate rows (see \code{\link{reshape}}).

  The \emph{long format} is called as such, because more sampling results in a
  lengthening of the data (i.e., the addition of rows) without any widening
  (i.e., the addition of columns). In contrast, it is common in community
  ecology to store abundance data as spreadsheets with sites as
  rows and species as columns. Such a data storage format is often called
  \emph{wide format}, because more sampling may result additionally in a
  widening of the data (e.g., more columns are required as further sampling
  reveals an increasing number of species).

  Function \code{stackSpecies} converts a matrix like input object, or file,
  into a representation in \emph{long format}, where each observation on a
  species consists of one row in a resulting data frame that has as many rows as
  there are single species occurrences in the data set, including layer
  replicates (see \code{\linkS4class{Species}}).
  
  In a similar fashion, function \code{stackSites} can be used to derive the
  \emph{long format} from a plot by variable matrix
  (see \code{\linkS4class{Sites}}).
  
  Function \code{stackCoordinates} collects spatial coordinates from spatial
  data stored on disk by calls to \code{\link[rgdal]{readOGR}}. It then reshapes
  two or three columns (in case of 3D features) of spatial coordinates along
  with a plot identifier to the format defined by class \verb{Sites}.
  
  Results of \code{stackSites} and \code{stackCoordinates} can be combined using
  \code{bind} to associate the sites with spatial information as intended for
  \verb{Vegsoup} objects (see \sQuote{Examples}).
}

\usage{
stackSpecies(x, file, sep = ";", dec = ",", schema = c("abbr", "layer"),
             discard = c("taxon", "comment"), absences, zeros = FALSE, verbose = FALSE)

stackSites(x, file, sep = ";", dec = ",", schema = "plot", zeros = FALSE, verbose = FALSE)

stackCoordinates(dsn, layer, schema, round = TRUE, verbose = TRUE, ...)
}

\arguments{
# Arguments to all functions
  \item{ x }{
     An object converted to a data frame.
  }
  \item{ file }{
     character. Path to a csv-file
  }
  \item{ sep, dec }{
    character. See \code{\link{read.csv}}.
  }
  \item{ schema }{
     character. Names of columns to be searched for. For \code{stackSites()} \code{schema = "rownames"} is available. See \sQuote{Details}.
  }
  \item{ zeros }{
     logical. If \code{TRUE}, drop leading zeros in plot names.
  }
  \item{ verbose }{
     Prints some diagnostic messages.
  }
# Arguments for function stackSpecies
  \item{ absences }{
    Character used to code absences. Can be missing, see \sQuote{Details}.
  }
  \item{ discard }{
     Names of additinal columns to be discarded.
  }
# Arguments for function stackCoordinates
  \item{ dsn }{
    Data source name. See \code{\link[rgdal]{readOGR}}.
  }
  \item{ layer }{
    Layer name. See \code{readOGR}.
  }
  \item{ round }{
    Round decimals to given precision. See \sQuote{Details}.
  }
  \item{ \dots }{
    additional arguments passed to \code{readOGR}. For \code{bind} method, objects of class
    \code{Species} respectively \code{Sites}.
  }
}

\details{
  For \code{stackSpecies} the supported data frame, either read from file, or
  passed as \code{R} object must hold columns corresponding to argument schema.
  The default schema is 'abbr', 'layer' (cf. class
  \code{\linkS4class{Species}}). The first element of schema has to name a
  column holding valid strings for 'abbr'. The second element defines the column
  that has layer memberships. Any further elements (e.g 'taxon') are only
  mandatory if \code{x} has additional columns except schema and the taxa block
  (the plain matrix of only species cover values, see below). These columns may
  also be covered by argument \code{discard} (possibly by its defaults).
  Any further fields can be used to store taxonomic or other relevant
  information (accuracy of determination) on a particular observation, again,
  this is not mandatory. The species matrix is assumed to have species in rows
  and plots in columns. Plot names are derived form \code{names(x)}. If argument
  file is supplied, data is imported using
  \code{read.csv(..., colClasses = "character")}. This matters if species
  abundances were recorded on a non-ordinal scale (e.g. integers) and shall
  later be treated as continuos. The \code{\linkS4class{Coverscale}}) class
  controls these conversions.

  Internal checking of type conversion is performed using
  \code{\link{type.convert}} to test for appropriate mode of the species
  abundances vector; setting \code{verbose = TRUE} will report results of this
  test. If argument \code{absences} is missing an almost save guess
  (matrix fill lower than 50\%) !is obtained
  from the data and derived as the most frequent value in the species data
  matrix, these are usually species absences.

  For \code{stackSites} there are no specific requirements about the form of the
  object, except that one column must be specified which will then be used as
  the plot identifier. The default is to search for a column named
  \code{'plot'}, but an other column name can be specified by argument
  \code{schema}. 

  For \code{stackCoordinates} the rgdal package needs to be loaded first.
  The first element of argument schema defines the name of the column in the OGR
  data source that identifies unique plots. This argument is mandatory and must
  match a column name in the spatial attributes data. A second element of
  argument \code{schema} can specify the name of the column in the OGR data
  source that identifies altitude measurements. If absent and the geometry is
  3D the function reads the third (z) dimension of the data source. If
  \code{verbose} the function will print what is returned by
  \code{\link[rgdal]{ogrInfo}}. Irrespective of the geometry type the function
  will use the \code{\link[sp:coordinates]{coordinates}} method in package
  \pkg{sp} to obtain a pair of geographic coordinates for each plot.
  In case of irregular polygons the centroid so obtained must not be accurate.
  If the coordinate reference system is not equivalent to
  \code{CRS("+init=epsg:4326")} it will be transformed.
}
\note{
  Functionality related to package \pkg{rgdal} is only available if this package
  is installed and loaded.
}
\value{
\code{stackSpecies} returns an object of class \code{\linkS4class{Species}}
\code{stackSites} and \code{stackCoordinates} returns an object of class \code{\linkS4class{Sites}}.
}

\author{ Roland Kaiser }

\examples{
# species matrix from the dune data set
require(vegan)
data(dune)

# create data.frame according to the defaults of the schema argument in stackSpecies
# there are two moss species in the dune data set, we assign them to a moss layer
x <- data.frame(abbr = names(dune),
                layer = c(rep("hl", 8), "ml",
                          rep("hl", 6), "ml",
                          rep("hl", 14)),
                t(dune))

# promote to class 'Species'
X <- stackSpecies(x, verbose = TRUE)

# stack corresponding sites data
data(dune.env)
x <- data.frame(plot = row.names(dune.env), dune.env)
Y <- stackSites(x, verbose = TRUE) # default schema = "plot" 

# the same using schema = "rownames"
Y <- stackSites(dune.env, schema = "rownames")

## create a file of dummy coordiantes in a common format for spatial data
# some random points
Y2 <- data.frame(x = rnorm(20), y = rnorm(20),
                  plot = unique(X$plot))

# promate to class 'SpatialPointsDataFrame'
coordinates(Y2) <- ~x+y
proj4string(Y2) <- CRS("+init=epsg:4326") # WGS84

# for demonstartion save to ESRI Shapefile
require(rgdal)
dsn <- tempfile()
layer <- "foo"
writeOGR(Y2, dsn, layer, driver = "ESRI Shapefile")

# read this file and promote to class 'Sites'
Y2 <- stackCoordinates(dsn, layer, schema = "plot")

# bind with remaining environmental data
Y <- bind(Y, Y2)
variables(Y)
}

\seealso{
\code{\link{bind}},
\code{\link{Species-class}},
\code{\link{Sites-class}},
\code{\link{reshapeSpecies}},
\code{\link{elevation}},
\code{\link{tile2latlng}},
\code{\link{read.csv}},
\code{\link[rgdal]{readOGR}} in package \pkg{rgdal}
}

\keyword{ import }

