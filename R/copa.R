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

## A function to create a matrix of permuted group assignments

perm.mat <- function(B, ids){
	mat <- matrix(NA, nrow=B, ncol=length(ids))
	for(i in 1:B){
		mat[i,] <- sample(ids, length(ids))
	}
	invisible(return(mat))
}


copa <- function(object, cl, cutoff = 5, max.overlap = 0, norm.count = 0, pct = 0.95){
  
  mat <- copaFilter(object, cl, cutoff, norm.count, pct)
  
  ng <- dim(mat)[1]
  ncomp <- choose(ng, 2)
  if(dim(mat)[1] > 1000){
    doit <-  getans(paste("\nYou are attempting to run copa on a large number of genes(",
                          dim(mat)[1], ").\n Which will result in ", ncomp, " pairwise comparisons.\n",
                          "This may take a long time and/or require lots of memory.\n",
                          "An alternative would be to increase the value of pct higher than",
                          pct, ".\n", "Do you want to continue?", sep=""))
    
    if(doit == "n"){
      return(cat("copa aborted at user request\n"))
    }
  }
  
  ## Get gene pairs that have mutually exclusive outliers
  outlier <- mat[,cl == 2] > cutoff 
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
                        mat = mat, cl = cl, cutoff = cutoff,
                        max.overlap = max.overlap, norm.count = norm.count,
                        pct = pct), class = "copa"))
}

plot.copa <- function(copa, idx, lib = NULL, sort = TRUE, col = NULL){

  ## If affy chip, use available library
  if(!is.null(lib))
    require(lib, character.only = TRUE) || stop(paste("The package", lib, "is not installed!\n"))
  
  layout(matrix(1:2, nc=1))
  ## convert ordered pair matrix into vector
  if(is.matrix(copa$ord.pr[idx,]))
    gn.idx <- as.vector(apply(copa$ord.prs[idx,], 1, c))
  else
    gn.idx <- copa$ord.prs[idx,]
  grps <- unique(sort(copa$cl))
  for(i in seq(along = gn.idx)){
    if(!is.null(lib))
      symb <- get(row.names(copa$mat)[gn.idx[i]],
                  get(paste(lib, "SYMBOL", sep="")))
    else
      symb <- row.names(copa$mat)[gn.idx[i]]
    if(sort){
      tmp.lst <- vector("list", length(grps))
      for(j in seq(along = grps)){
        ord <- order(unlist(copa$mat[gn.idx[i], copa$cl == grps[j]]))
        tmp.lst[[j]] <- copa$mat[gn.idx[i], copa$cl == grps[j]][ord]
      }
      barplot(unlist(tmp.lst), main = symb, col = rep(col, table(copa$cl)))
    }else{
      tmp.lst <- vector("list", length(grps))
      for(j in seq(along = grps))
        tmp.lst[[j]] <- copa$mat[gn.idx[i], copa$cl == grps[j]]
        barplot(unlist(tmp.lst), main = symb, col = rep(col, table(copa$cl)))
    }
  }
}



copaPerm <- function(object, copa, outlier.num, B = 100, pval = FALSE){
 
  perm <- perm.mat(B, copa$cl)
  prmvals <- vector("list", B)
  for(i in 1:B){
    mat <- copaFilter(object, copa$cl, copa$cutoff, copa$norm.count, copa$pct)

    ## Get gene pairs that have mutually exclusive outliers
    outlier <-mat[,perm[i,] == 2] > copa$cutoff 
    outpairs <- crossprod(t(outlier))

    ## Get row sums and pairwise row sums
    r.sum <- rowSums(outlier)
    pr.sums <- pSum(r.sum)
    
    ## find pairs where outliers are mutually exclusive and row sums are 'large'
    ## use lower.tri to remove diagonal and duplicated comparisons
    outpairs <- outpairs[lower.tri(outpairs)]
    pr.sums <- pr.sums[lower.tri(pr.sums)]
    index <- outpairs <= copa$max.overlap
    prmvals[[i]] <- pr.sums[index]
  } 
  out <- sapply(prmvals, function(x) sum(x >= outlier.num))
  if(pval){
    p.value <- sum(out <= outlier.num)/B
    fdr <- mean(out)
    return(list(out = out, p.value = p.value, fdr = fdr))
  }else{
    return(out)
  }
}

copaFilter <- function(object, cl, cutoff, norm.count, pct){
  
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

  ## Median center and scale using MAD
  med <- apply(mat, 1, median)
  MAD <- apply(mat, 1, mad)
  mat <- sweep(mat, 1, med, "-")
  mat <- sweep(mat, 1, MAD, "/")

  outlier <- mat[,cl == 2] > cutoff
  badnorm <- mat[,cl == 1] > cutoff

  num.out <- rowSums(outlier)
  num.bad <- rowSums(badnorm)
  
  num.out <- num.out[num.bad <= norm.count]
  mat <- mat[num.bad <= norm.count,]
  index <- num.out >= quantile(num.out, probs = pct)
  mat <- mat[index,]
  mat
}
  
