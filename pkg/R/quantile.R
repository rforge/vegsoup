#	summary statistics
#	replaces fivenum (Tukey Five-Number Summary)
#if (!isGeneric("quantile")) {
setGeneric("quantile",
	function (x, ...)
		standardGeneric("quantile")
)
#}
.quantile.Vegsoup <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7, coverscale = FALSE, ...) {
	if (coverscale & !is.null(decostand(x))) {
		message("disregard decostand method for calculations of quantiles")
	 	decostand(x) <- NULL
	}
	if (!inherits(x, "VegsoupPartition")) {
		#	to keep the function below unchanged
		x <- VegsoupPartition(x, clustering = rep(1, nrow(x)))
	}
	xx <- as.numeric(x)
	xx[ xx == 0 ] <- NA
	tmp <- aggregate(xx,
		by = list(partitioning(x)),
		FUN = function (x) stats::quantile(x, probs = probs, na.rm = TRUE), simplify = FALSE) # , ...
	part <- tmp[ , 1 ]
	tmp <- tmp[ , -1 ]
	res <- array(0, dim = c(dim(tmp)[ 2 ], dim(tmp)[ 1 ], length(probs)),
		dimnames = list(names(tmp), part, probs))
	for (i in seq(along = probs)) {
		for (j in 1:nrow(res)) {
			#	j = 1; i = 1
			res[ j, , i ] <- sapply(tmp[ , j ], function (x) x[ i ])
		}
	}
	#	groome names
	dimnames(res)[[3]] <- paste0("q", dimnames(res)[[3]])
	if (coverscale & is.ordinal(x)) {
		#	message("backconvert to coverscale: ", coverscale(x)@name)
		for (i in 1:dim(res)[3]) {
			tmp <- res[, , i]
			if (mode(tmp) != "numeric") mode(tmp) <- "numeric"
			vals <- as.character(cut(tmp,
				breaks = c(-1, coverscale(x)@lims, 100),
				labels = c(".", coverscale(x)@codes), right = FALSE))
			tmp.i <- tmp
			mode(tmp.i) <- "character"
			tmp.i[] <- vals
			res[, , i] <- tmp.i
		}
	} else {
			#	message("did not backconvert to coverscale: ", coverscale(x)@name)	
		if (coverscale & is.continuous(x)) {
			message("coverscale is not ordinal")
		}
	}
	res
}

setMethod("quantile",
	signature(x = "VegsoupPartition"),
	function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7, coverscale = FALSE, ...) {
		.quantile.Vegsoup(x, probs = probs, na.rm = na.rm, names = names, type = type, coverscale = coverscale, ...)	
	}
)

setMethod("quantile",
	signature(x = "Vegsoup"),
	function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7, coverscale = FALSE, ...) {
		.quantile.Vegsoup(x, probs = probs, na.rm = na.rm, names = names, type = type, coverscale = coverscale, ...)	
	}
)