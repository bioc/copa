###############################################
##
## copa package copyright 2006 James W. MacDonald
##
## copa - a package to implement the cancer outlier profile analysis
##        of Tomlins, et. al. Science. 2005 Oct 28;310(5748):644-8.
##
## November 2005 - original implementation as separate functions
##
## January 2006 - converted to package form and generalized to any
##                microarray platform
##
#################################################

## Function to get all pairwise sums and return in a matrix

pSum <- function(a){
  if(any(is.na(a)))
    stop(paste("The centered and scaled matrix", mat, "cannot contain any NA values!\n",
               "This is often due to genes with identical values for all samples.\n",
               "Please filter these data first.\n"))
  out <- .C("pairwise_sum", as.double(a), as.integer(length(a)),
            as.double(vector("double", length(a)^2)))[[3]]
  matrix(out, nc = length(a), nr = length(a))
}

## Simple function to query end user
getans <- function(msg, allowed = c("y","n")){
  repeat{
    out <- paste("[", paste(allowed, collapse = "/"),
                 "]")
    outmsg <- paste(msg, out)
    cat(outmsg)
    answer <- readLines(n = 1)
    if(answer %in% allowed)
      break
    else cat(paste(answer, "is not a valid response, try again.\n"))
  }
  answer
}


copa <- function(object, cutoff = 5, max.overlap = 0){

  mat <- NULL
  if(is(object, "exprSet")){
    mat <- exprs(object)
  }else{
    if(is(object, "PLMset")){
      mat <- coefs.probe(object)
    }else{
      if(is(object, "MAList")){
        mat <- object$M
      }else{
        if(is(object, "marrayNorm")){
          mat <- object@maM
        }
      }
    }
  }
  if(is.null(mat))
    mat <- as.matrix(object)
  
  ng <- dim(mat)[1]
  ncomp <- ng * (ng - 1)/2
  if(dim(mat)[1] > 1000){
   doit <-  getans(paste("\nYou are attempting to run copa on a large number of genes(",
                 dim(mat)[1], ").\n Which will result in ", ncomp, " pairwise comparisons.\n",
                 "This may take a long time and/or require lots of memory.\n",
                 "Do you want to continue?", sep=""))
  
   if(doit == "n"){
     return(cat("copa aborted at user request\n"))
   }
 }
  
  ## Median center and scale using MAD
  med <- apply(mat, 1, median)
  MAD <- apply(mat, 1, mad)
  mat <- sweep(mat, 1, med, "-")
  mat <- sweep(mat, 1, MAD, "/")

  ## Get gene pairs that have mutually exclusive outliers
  outlier <- mat > cutoff
  outpairs <- crossprod(t(outlier))

  ## Get row sums and pairwise row sums
  r.sum <- rowSums(outlier)
  pr.sums <- pSum(r.sum)

  ## find pairs where outliers are mutually exclusive and row sums are 'large'
  ## use lower.tri to remove diagonal and duplicated comparisons
  outpairs <- outpairs[lower.tri(outpairs)]
  pr.sums <- pr.sums[lower.tri(pr.sums)]
  index <- outpairs <= max.overlap
  ord <- order(pr.sums[index], decreasing = TRUE)
  mat1 <- matrix(1:dim(mat)[1], nc = dim(mat)[1], nr = dim(mat)[1])
  mat2 <- t(mat1)           
  pr.idx <- cbind(mat1[lower.tri(mat1)], mat2[lower.tri(mat2)])
  ord.prs <- pr.idx[index,][ord,]
  return(structure(list(ord.prs = ord.prs, pr.sums = pr.sums[index][ord],
                        mat = mat), class = "copa"))
}

plot.copa <- function(copa, idx, cl, lib = NULL, sort = TRUE, col = NULL){

  ## If affy chip, use available library
  if(!is.null(lib))
    require(lib, character.only = TRUE) || stop(paste("The package", lib, "is not installed!\n"))
  
  layout(matrix(1:2, nc=1))
  ## convert matrix into vector ordered by row
  if(is.matrix(copa$ord.pr[idx,]))
    gn.idx <- as.vector(apply(copa$ord.prs[idx,], 1, c))
  else
    gn.idx <- copa$ord.prs[idx,]
  if(is.factor(cl)) grps <- levels(cl) else grps <- unique(cl)
  for(i in seq(along = gn.idx)){
    if(!is.null(lib))
      symb <- get(row.names(copa$mat)[gn.idx[i]],
                  get(paste(lib, "SYMBOL", sep="")))
    else
      symb <- row.names(copa$mat)[gn.idx[i]]
    if(sort){
      tmp.lst <- vector("list", length(grps))
      for(j in seq(along = grps)){
        ord <- order(unlist(copa$mat[gn.idx[i], cl == grps[j]]))
        tmp.lst[[j]] <- copa$mat[gn.idx[i], cl == grps[j]][ord]
      }
      barplot(unlist(tmp.lst), main = symb, col = rep(col, table(cl)))
    }else{
      tmp.lst <- vector("list", length(grps))
      for(j in seq(along = grps))
        tmp.lst[[j]] <- copa$mat[gn.idx[i], cl == grps[j]]
        barplot(unlist(tmp.lst), main = symb, col = rep(col, table(cl)))
    }
  }
}
