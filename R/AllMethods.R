## S4 methods for the copaFilter function

setGeneric("copaFilter", function(object, ...) standardGeneric("copaFilter"))

setMethod("copaFilter", "matrix", function(object, cl, cutoff, norm.count, pct){
  do.copaFilter(object, cl, cutoff, norm.count, pct)
})

setMethod("copaFilter", "data.frame", function(object, cl, cutoff, norm.count, pct){
  object <- as.matrix(object)
  do.copaFilter(object, cl, cutoff, norm.count, pct)
})

setMethod("copaFilter", "ExpressionSet", function(object, cl, cutoff, norm.count, pct){
  object <- exprs(object)
  do.copaFilter(object, cl, cutoff, norm.count, pct)
})

