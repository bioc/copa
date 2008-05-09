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
            as.double(vector("double", length(a)^2)), PACKAGE = "copa")[[3]]
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
  
  if(any(is.na(match(unique(cl), 1:2)))) stop("The 'cl' vector can only contain 1's and 2's",
                                              call. = FALSE)
  
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

  ## Get quantile differences for ordering ties
  q.norm <- apply(mat[,cl == 1], 1, quantile, probs = 0.75)
  q.tum <- apply(mat[,cl == 2], 1, quantile, probs = 0.75)
  q.diff <- q.tum - q.norm
  q.sums <- pSum(q.diff)

  ## find pairs where outliers are mutually exclusive and row sums are 'large'
  ## use lower.tri to remove diagonal and duplicated comparisons
  outpairs <- outpairs[lower.tri(outpairs)]
  pr.sums <- pr.sums[lower.tri(pr.sums)]
  q.sums <- q.sums[lower.tri(q.sums)]
  index <- outpairs <= max.overlap
  ord <- order(pr.sums[index], q.sums[index], decreasing = TRUE)
  mat1 <- matrix(1:dim(mat)[1], nc = dim(mat)[1], nr = dim(mat)[1])
  mat2 <- t(mat1)           
  pr.idx <- cbind(mat1[lower.tri(mat1)], mat2[lower.tri(mat2)])
  ord.prs <- pr.idx[index,][ord,]
  return(structure(list(ord.prs = ord.prs, pr.sums = pr.sums[index][ord],
                        mat = mat, cl = cl, cutoff = cutoff,
                        max.overlap = max.overlap, norm.count = norm.count,
                        pct = pct), class = "copa"))
}

plotCopa <- function(copa, idx, lib = NULL, sort = TRUE, col = NULL,
                     legend = NULL){

  ## If affy chip, use available library
  if(!is.null(lib))
    require(lib, character.only = TRUE) || stop(paste("The package", lib, "is not installed!\n"))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
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
      barplot(unlist(tmp.lst), main = symb, col = rep(col, table(copa$cl)),
              xaxt = "n")
      if(!is.null(legend))
        if(i%%2 == 1)
          legend("topleft", inset = 0.01, legend = legend, fill = col)
    }else{
      tmp.lst <- vector("list", length(grps))
      for(j in seq(along = grps))
        tmp.lst[[j]] <- copa$mat[gn.idx[i], copa$cl == grps[j]]
        barplot(unlist(tmp.lst), main = symb, col = rep(col, table(copa$cl)),
                xaxt = "n")
    }
  }
}



copaPerm <- function(object, copa, outlier.num, gene.pairs, B = 100,
                     pval = FALSE, verbose = TRUE){
 
  perm <- perm.mat(B, copa$cl)
  prmvals <- vector("list", B)
  if(verbose)  cat("Counting permutations...\n")
  mat <- copaFilter(object, copa$cl, copa$cutoff, copa$norm.count, copa$pct)
  for(i in 1:B){
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
    if(verbose)
      if(i %% 100 == 0) cat(paste(i, "\n"))
  } 
  out <- sapply(prmvals, function(x) sum(x >= outlier.num))
  if(pval){
    p.value <- sum(out >= gene.pairs)/B
    fdr <- mean(out)/gene.pairs * 100
    return(list(out = out, p.value = p.value, fdr = fdr))
  }else{
    return(out)
  }
}

do.copaFilter <- function(object, cl, cutoff, norm.count, pct){
  
  ## Median center and scale using MAD
  med <- apply(object, 1, median)
  MAD <- apply(object, 1, mad)
  ## Protect against samples with all equal values
  if(any(MAD == 0))
    MAD <- ifelse(MAD == 0, 0.001, MAD)
  object <- sweep(object, 1, med, "-")
  object <- sweep(object, 1, MAD, "/")

  outlier <- object[,cl == 2] > cutoff
  badnorm <- object[,cl == 1] > cutoff

  num.out <- rowSums(outlier)
  num.bad <- rowSums(badnorm)
  
  num.out <- num.out[num.bad <= norm.count]
  object <- object[num.bad <= norm.count,]
  if(quantile(num.out, probs = pct) == 0)
    stop(paste("The ", pct, " quantile is 0.\n",
               "You need to increase this value for",
               " copa to work correctly."))
  index <- num.out >= quantile(num.out, probs = pct)
  object <- object[index,]
  object
}
  
tableCopa <- function(copa){
  rev(table(copa$pr.sums))
}

scatterPlotCopa <- function(copa, idx, lib = NULL){
  if (!is.null(lib)) 
        require(lib, character.only = TRUE) || stop(paste("The package", 
            lib, "is not installed!\n"), call. = FALSE)
  for(i in seq(along = idx)){
    gn.idx <- copa$ord.pr[idx[i],,drop = FALSE]
    if(!is.null(lib)){
      xlab <- get(row.names(copa$mat)[gn.idx[i,1]],
                  get(paste(lib, "SYMBOL", sep = "")))
      ylab <- get(row.names(copa$mat)[gn.idx[i,2]],
                  get(paste(lib, "SYMBOL", sep = "")))
    }else{
      xlab <- row.names(copa$mat)[gn.idx[i,1]]
      ylab <- row.names(copa$mat)[gn.idx[i,2]]
    }
    plot(copa$mat[gn.idx[i,1],], copa$mat[gn.idx[i,2],], xlab = xlab,
         ylab = ylab, col = copa$cl)
  }
}
  
summaryCopa <- function(copa, pairnum, lib = NULL){
    if(pairnum > max(copa$pr.sums))
        stop("There aren't any genes with that many pairs.\n", call. = FALSE)
    if(!is.null(lib))
        require(lib, character.only = TRUE) || stop(paste("The package",
                     lib, "is.not installed!\n"), call. = FALSE)
    idx <- copa$pr.sums >= pairnum
    prbId1 <- row.names(copa$mat)[copa$ord.prs[idx,1]]
    prbId2 <- row.names(copa$mat)[copa$ord.prs[idx,2]]
    if(is.null(lib)){
        out <- data.frame("Number of pairs" = copa$pr.sums[idx],
                          "Probe ID 1" = prbId1,
                          "Probe ID 2" = prbId2)
    }else{
        gn1 <- sapply(mget(prbId1, get(paste(lib, "SYMBOL", sep = ""))),
                           function(x) x[1])
        gn2 <- sapply(mget(prbId2, get(paste(lib, "SYMBOL", sep = ""))),
                           function(x) x[1])
        out <- data.frame("Number of pairs" = copa$pr.sums[idx],
                          "Probe ID 1" = prbId1,
                          "Symbol 1" = gn1,
                          "Probe ID 2" = prbId2,
                          "Symbol 2" = gn2, row.names = NULL)
    }
    out
}
