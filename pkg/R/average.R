setGeneric("average",
	function (obj, ...)
		standardGeneric("average")
)

setMethod("average",
	signature(obj = "VegsoupPartition"),
	function (obj, ...) {
		cl <- match.call()
		if (any(names(cl) == "mode")) {
			if (cl$mode == "R") {
				res <- t(aggregate(as.matrix(obj, ...),
					by = list(splitAbbr(obj)$layer),
					FUN = sum)[, -1])
				res <- res / contingency(obj, ...)
				#	division by zero
				res[is.na(res)] <- 0
				colnames(res) <- layers(obj)
			} else {
				stop("omit mode argument for standard behaviour")
			}
		} else {
			res <- t(aggregate(as.matrix(obj, ...),
				by = list(partitioning(obj)), FUN = sum))[-1, ]
			res <- res / contingency(obj)
			#	division by zero
			res[is.na(res)] <- 0
			colnames(res) <- unique(partitioning(obj))
			rownames(res) <- colnames(obj)
		}
		return(res)
	}
)