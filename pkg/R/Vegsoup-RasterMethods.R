#if (!isGeneric("extract")) {
	setGeneric("extract", function (x, y, ...)
		standardGeneric("extract"))
#}

setMethod("extract",
	signature(x = "Raster", y = "Vegsoup"), 
	function(x, y, ...){
		require(raster) 
		raster::extract(x, as(y, "SpatialPoints"), ...)
	}
)