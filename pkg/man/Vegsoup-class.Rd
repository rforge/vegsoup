\encoding{UTF-8}

\name{Vegsoup-class}
\Rdversion{1.1}

\docType{class}

\alias{Vegsoup-class}

\title{ Class Vegsoup }

\description{
  This is the foundation class of all \emph{Vegsoup*} objects.
}

\section{ Objects from the Class }{

Objects can be created by calls of the form \code{new("Vegsoup", ...)}.  
The user is encouraged to use the creator function \code{\link{Vegsoup}}.  

}
\section{ Slots }{
  \describe{
	\item{\code{species}:}{
	  Object of class \code{\linkS4class{Species}},
	  the species data in \emph{long} format.
	}
	\item{\code{sites}:}{
	  Object of class \code{\linkS4class{data.frame}},
	  the sites data in \emph{wide} format, developed from
	  class \code{\linkS4class{Sites}} in \emph{long} format.
	}
	\item{\code{taxonomy}:}{
	  Object of class \code{\linkS4class{Taxonomy}},
	  the taxonomic reference list.
	}
	\item{\code{coverscale}:}{
	  Object of class \code{\linkS4class{Coverscale}},
	  the description of the cover scale. See
	  \code{\linkS4class{Coverscale}} for details.
	}
	\item{\code{decostand}:}{
	  \code{character}, % 
	  the method for \code{\link[vegsoup]{decostand}} in \pkg{vegan}. 
	  Takes effect if \code{\link{as.matrix}} or \code{\link{as.dist}}
	  is called. The \emph{slot 'decostand'} is in fact of class
	  \code{decostand}, this class is not exposed to the user
	  and hence not documented.
	}
	\item{\code{dist}:}{
	  \code{character},
	  giving a method for \code{\link[vegsoup]{vegdist}}.
	  Takes effect in calls or dispatch of \code{\link{as.dist}}.
	}
	\item{\code{layers}:}{
	  \code{character},
	  containing the (vegetation) layers the object is structured by.
	  Order matters.
	}
	\item{\code{group}:}{
	  \code{integer},
	  a vector of a-priori group assignments.
	}
	\item{\code{sp.points}:}{
	  Object of class \code{\linkS4class{SpatialPointsDataFrame}}
	  defined in package \pkg{sp}. See also \code{\link{SpatialPoints}}
	}
	\item{\code{sp.polygons}:}{
	  Object of class \code{\linkS4class{SpatialPolygonsDataFrame}}.
	  See also \code{\link{SpatialPolygons}}.
	}
  }
}

\section{ Methods }{

  Classes \code{\linkS4class{VegsoupOptimstride}},
  \code{\linkS4class{VegsoupPartition}} and \code{"VegsoupPartitionFidelity"}
  all contain and directly extend class \verb{Vegsoup}. The methods
  presented in the following sections apply to all these classes.

  \describe{
%% 2
  \bold{ 1) Dimensions of the object (species matrix).}
	\item{dim}{
	  \code{signature(x = "Vegsoup")}:
	  returns matrix dimensions (see above) giving the number of plots and
	  (pseudo-)species.
	}
	\item{ncol}{
	  \code{signature(x = "Vegsoup")}:
	  the number of (pseudo-)species.
	}
	\item{nrow}{
	  \code{signature(x = "Vegsoup")}:
	  the number of plots.
	}
	\item{ncell}{
	  \code{signature(x = "Vegsoup")}:
	  the number of cells in the species matrix.
	}
	\item{MatrixFill}{
	  \code{signature(obj = "Vegsoup")}:
	  the percentage of presences in the species matrix.
	}

%% 2
  \bold{ 2) Layer structure.}

	\item{layer}{
	  \code{signature(obj = "Vegsoup")}:
	  returns the layer for each species observation.
	}
	\item{layers}{
	  \code{signature(obj = "Vegsoup")}:
	  modifies an object's layer structure by suppling a layer combination
	  rule as a \code{"character"} vector \code{'collapse'}.
	}
	\item{layers<-}{
	  \code{signature(obj = "Vegsoup")}:
	  modifies the order of layers. This has effects on \code{as.matrix} and
	  similar.
	}

%% 3
  \bold{ 3) Cover scale.}
	\item{coverscale}{
	  \code{signature(obj = "Vegsoup")}:
	  returns the cover scale.
	}
	\item{coverscale<-}{
	  \code{signature(obj = "Vegsoup", value = "character")}:
	  sets the coverscale by suppling a name as recognized by the class
	  creator function \code{\link{Coverscale}}.
	  \code{signature(obj = "Vegsoup", value = "Coverscale")}:
	  sets the coverscale by suppling a \code{"Coverscale"} object.
	}
	\item{BraunBlanquetReduce}{
	  \code{signature(obj = "Vegsoup")}:
		 reduce 9 point Braun-Blanquet scale to 7 point.
	}

%% 4
  \bold{ 4) Standardization methods.}
	\item{decostand}{
	   \code{signature(obj = "Vegsoup")}:
	   returns the standardization method active for an object.
	}
	\item{decostand<-}{
	  \code{signature(obj = "Vegsoup", value = "character")}:
	   sets the standardization method active for an object. Note,
	   values of \emph{slot 'decostand'} have effect in \code{as.matrix} and
	   similar matrix methods.
	  \code{signature(obj = "Vegsoup", value = "NULL")}:
	   reset/remove the standardization for an object by setting the value to
	   \code{NULL}.
	}

%% 5
  \bold{ 5) Dissimilarity/distance index.}
	\item{vegdist}{
	  \code{signature(x = "Vegsoup")}:
	  returns the dissimilarity/distance index set for an object.
	}
	\item{vegdist<-}{
	  \code{signature(x = "Vegsoup", value = "character")}:
	  set a distance measure for an object.
	}
	\item{distconnected}{
	  \code{signature(x = "Vegsoup")}:
	  connectedness of Dissimilarities.
	}
%% 6
  \bold{ 6) Retrieve indices to species matrix.}
	\item{indices}{
	  \code{signature(x = "Vegsoup")}:
	  returns the locations (row and column pointers) and values of all
	  nonzero entries in a species matrix (see \code{as.matrix} method and
	  alike below).
	}
	\item{as.vector}{
	  \code{signature(x = "Vegsoup", mode = "missing")}:
	  is a shortcut to \code{as.vector(as.matrix(x, typeof = mode), mode)}.
	  Note, here the \code{'mode'} argument has the meaning of the base
	  \code{R} function and specifies the storage mode of the returned vector.
	 }
   Together, \code{indices} and \code{as.vector} can be used to create a
   species matrix. See \code{\link{as.matrix}}.

%% 7
  \bold{ 7) Retrieve species matrix.}
	In fact there are several ways of how to obtain the species matrix from an
	\code{Vegsoup*} object and all are based on the logic of object coercion.
	There are methods focused on a specific storage mode, as well, as a
	generic \code{as.matrix} coercion method. The following section outlines
	even more general index methods.
	\item{as.character}{
	  \code{signature(x = "Vegsoup")}:
	  returns a matrix of mode "character" with species abundances in the data
	  set's original scale.
	}
	\item{as.logical}{
	  \code{signature(x = "Vegsoup")}:
	  returns a matrix of mode "numeric" with zeros coding for species
	  absences and ones for presences.
	}
	\item{as.numeric}{
	  \code{signature(x = "Vegsoup")}:
	  returns a matrix of mode "numeric" with species abundances recoded
	  according to \code{coverscale(obj)}.
	}
	\item{as.matrix}{
	  \code{signature(x = "Vegsoup")}:
		retrieves a species matrices like those outlined above by setting the
		\code{typeof} argument specifying the mode of the returned matrix; can
		be one of \code{"character"}, \code{"numeric"} or \code{"logical"}.
		Use the additional argument \code{mode} with all four methods to get
		matrix transpositions. See \code{\link{as.matrix}}.
	}

%% 8
  \bold{ 8) Retrieve dissimilarity/distance matrix.}
	\item{as.dist}{
	  \code{signature(m = "Vegsoup")}:
	  compute a distance matrix as specified by \code{vegdist(obj)}.
	  See \code{\link{as.dist}}
	}
	\item{getDistconnected}{
	  \code{signature(obj = "Vegsoup")}:
	  tests if dissimilarities/distances are connected.
	}

%% 9
  \bold{ 9) Retrieve class slots.}
	Methods to access or replace the \code{'species'}, \code{'sites'}
	or \code{'taxonomy'} slots are:
	\item{species}{
	  \code{signature(obj = "Vegsoup")}:
	  returns slot, on object of class
	  \code{\linkS4class{Species}}.
	}
	\item{species<-}{
	  \code{signature(obj = "Vegsoup", value = "Species")}:
	  assigns object of class \verb{Species} to an existing
	  \code{Vegsoup*} object. See \code{\link{species<-}}.
	  \code{signature(obj = "Vegsoup", value = "SpeciesTaxonomy")}:
	  assigns object of class \code{"SpeciesTaxonomy"} to an existing
	  \code{Vegsoup*} object. Not implemented yet! See
	  \code{\linkS4class{SpeciesTaxonomy}}.
	}
	\item{sites}{
	  \code{signature(obj = "Vegsoup")}:
	  returns slot, an object of class
	  \code{"data.frame"}. Note, the sites slot of Vegsoup objects
	  is of class data.frame not \verb{Sites}.
	  %% equaling an object of class \code{\linkS4class{Sites}}.
	}
	\item{Sites<-}{
	  \code{signature(obj = "Vegsoup", value = "data.frame")}:
	  The value object for the replacement function needs to have rownames
	  matching \code{rownames(obj)}.
	}
	\item{Taxonomy}{
	  \code{signature(obj = "Vegsoup")}:
	  returns slot, an object of class
	  \code{"data.frame"}, equaling an object of class
	  \code{\linkS4class{Taxonomy}}.
	}
	\item{Taxonomy<-}{
	  \code{signature(obj = "Vegsoup", value = "SpeciesTaxonomy")}:
	  assign object of class \code{"SpeciesTaxonomy"} to an existing
	  \code{Vegsoup*} object. Not implemented yet! See
	  \code{\linkS4class{SpeciesTaxonomy}}
	}
	Methods to access or replace the \code{'Spatial*'} slots are:
	\item{SpatialPointsVegsoup}{
	  \code{signature(obj = "Vegsoup")}:
	  returns object of class \code{\linkS4class{SpatialPointsDataFrame}}.
	}
	\item{SpatialPolygonsVegsoup}{
	  \code{signature(obj = "Vegsoup")}:
	  returns object of class \code{\linkS4class{SpatialPolygonsDataFrame}}.
	}
	\item{apriori}{
	  \code{signature(obj = "Vegsoup")}:
	  returns a vector of a-priori group assignments. If not defined an
	  \code{"integer"} vector of length equal to \code{nrow(obj)}.
	}

%% 10
  \bold{ 10) Names and species abbreviations.}
	\item{dimnames}{
	  \code{signature(x = "Vegsoup")}:
	  dimnames of the species matrix. See \code{\link{dimnames}}.
	}
	\item{rownames}{
	  \code{signature(x = "Vegsoup")}:
	  rownames of the species matrix, as well as of the sites data frame.
	}
	\item{colnames}{
	  \code{signature(x = "Vegsoup")}:
	  the column names of the species \code{"matrix"}.
	  See \code{\link{abbr.layer}}.
	}
	\item{names}{
	  \code{signature(x = "Vegsoup")}:
	  the column names of the sites \code{"data.frame"}.
	}
	\item{names<-}{
	  \code{signature(x = "Vegsoup", value = "character")}:
	}
	\item{abbr}{
	  \code{signature(obj = "Vegsoup")}:
	  abbreviations used to encode species.
	}
	\item{abbr.layer}{
	  \code{signature(obj = "Vegsoup")}:
	   abbreviations used to encode (pseudo-)species (layer replicates).
	}
	\item{splitAbbr}{
	  \code{signature(obj = "Vegsoup")}:
	  returns a \code{data.frame}. See \code{\link{splitAbbr}}.
	}

%% 11
  \bold{ 11) Extract and replace.}
	\item{[}{
	  \code{signature(x = "Vegsoup", i = "ANY", j = "ANY")}:
	  provides all common indexing operations allowed for class
	  \code{"matrix"}.
	}
	\item{[<-}{
	  \code{signature(x = "Vegsoup", i = "ANY", j = "missing")}:
	  allows assignments to \emph{slot 'sites'} . Note, this does not
	  affect the species matrix.
	}
	\item{$}{
	  \code{signature(x = "Vegsoup")}:
	  access columns of sites data by name.
	}
	\item{$<-}{
	  \code{signature(x = "Vegsoup")}:
	  assign columns of sites data by name.
	}

%% 12
  \bold{ 12) Row and column aggregation.}
	\item{colSums}{
	  \code{signature(x = "Vegsoup")}:
	  returns \code{colSums(as.logical(obj))}, hence does not take
	  species abundance/cover values into consideration.
	}
	\item{rowSums}{
	  \code{signature(x = "Vegsoup")}:
	  returns \code{rowSums(as.logical(obj))}, same as above but for rows.
	}
	\item{colMeans}{
	  \code{signature(x = "Vegsoup")}:
	  returns \code{colMeans(as.numeric(obj))} of species
	  abundance/cover values.
	}
	\item{rowMeans}{
	  \code{signature(x = "Vegsoup")}:
	  returns \code{rowMeans(as.numeric(obj))}.
	}
	\item{richness}{
	  \code{signature(obj = "Vegsoup")}:
	  Conceptually similar to the methods outlined above, but tailored to
	  respect the layer structure of an object in a special way.
	  See \code{\link{richness}} for details.
	}

%% 13
  \bold{ 13) Reorder and sample.}
	\item{seriation}{
	  \code{signature(obj = "Vegsoup")}:
	  modifies order as to reflect a vegetation table.
	}
	\item{sample}{
	  \code{signature(x = "Vegsoup")}:
	  random samples and permutations.
	}


%% 14
  \bold{ 14) Plotting methods.}
	\item{plot}{
	  \code{signature(x = "Vegsoup", y = "missing")}:
	  a simple descriptive plot of the data set.
	}
	\item{coldiss}{
	  \code{signature(obj = "Vegsoup")}:
	  plots of the dissimilarity matrix.
	}
	\item{map}{
	  \code{signature(obj = "Vegsoup")}:
	  plots a plain map.
	}

%% 15
  \bold{ 15) Spatial methods.}
	\item{coordinates}{
	  \code{signature(obj = "Vegsoup")}:
	  returns the a \code{"matrix"} of coordinates.
	  See \code{\link[sp:coordinates]{coordinates}}.
	}
	\item{elevation}{
	  \code{signature(obj = "Vegsoup")}:
	  queries elevation for plot coordinates from the SRTM3 data set.
	}
	\item{quicklook}{
	  \code{signature(obj = "Vegsoup")}:
	  open a map in browser window.
	}

%% 16
  \bold{ 16) Summary methods.}
	\item{summary}{
	  \code{signature(object = "Vegsoup")}:
	  prints an summary.
	}
	\item{head}{
	  \code{signature(x = "Vegsoup")}:
	  prints the top rows of the species matrix.
	}
	\item{tail}{
	  \code{signature(x = "Vegsoup")}:
	  prints the bottom rows of the species matrix.
	}
	\item{taxalist}{
	  \code{signature(obj = "Vegsoup")}:
	  returns a complete list of species with corresponding abbreviations.
	}

%% 17
  \bold{ 17) Coerce methods.}
	\item{coerce}{
	  \code{signature(from = "Vegsoup", to = "dist")}:
	  same as \code{as.dist(obj)} (see above).
	  \code{signature(from = "Vegsoup", to = "list")}:
	  returns objects all of class \code{"data.frame"} from slots
	  \code{'species'}, \code{'sites'} and \code{'taxonomy'}.
	  \code{signature(from = "Vegsoup", to = "matrix")}:
	  same as \code{as.matrix(obj)} (see above).
	  \code{signature(from = "Vegsoup", to = "data.frame")}:
	  returns objects all of class \code{"data.frame"}. See
	  \code{\link{as.data.frame}}. %% improve doc
	  \code{signature(from = "Vegsoup", to = "data.list")}:
	  returns objects all of class \code{"data.list"}. See
	  \code{\link{as.data.list}}. %% improve doc
	  \code{signature(from = "Vegsoup", to = "dsparseMatrix")}:
	  returns objects all of class \code{"dsparseMatrix"}. See
	  \code{\link{as.matrix}}. %% improve doc
	  \code{signature(from = "Vegsoup", to = "sparseMatrix")}:
	  returns objects all of class \code{"sparseMatrix"}. See
	  \code{\link{as.matrix}}. %% improve doc
	}

%% 18
  \bold{ 18) Combine objects.}
	\item{bind}{
	  \code{signature(... = "Vegsoup")}:
	  combine/merge objects.
	}

%% 19
  \bold{ 1) Output methods.}
	\item{KML}{
	  \code{signature(obj = "Vegsoup")}:
	  writes a KML file for visualization in Google Earth.
	}	
	\item{write.verbatim}{
	  \code{signature(obj = "Vegsoup")}:
	  writes the species matrix including header information to a plain text
	  file.
	}

%% 20
  \bold{ 1) Miscellaneous.}
	\item{Indpower}{
	  \code{signature(obj = "Vegsoup")}:
	  congruence between indicator and target species.
	}
	\item{compspec}{
	  \code{signature(obj = "Vegsoup")}:
	  compositional indicator species analysis.
	}
  }
}

\author{ Roland Kaiser }

\note{
Give notes!
}

\seealso{
\code{\linkS4class{Species}},
\code{\linkS4class{Sites}},
\code{\linkS4class{Taxonomy}},
\code{\linkS4class{Vegsoup}},
\code{\link{instantsoup}}
\code{\link{coenoflex}}
}
\examples{
showClass("Vegsoup")
}
\keyword{ classes }
\keyword{ methods }

