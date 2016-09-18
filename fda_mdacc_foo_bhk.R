########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 15, 2015
########################


## courtesy of Matthew McCall
celfileDateHour <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	if(length(h$ScanDate) > 0) {
	    h$ScanDate <- gsub(pattern="T", replacement=" ", x=h$ScanDate)
	    ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else { ddate <- rep(NA, 2)}
    names(ddate) <- c("day", "hour")
	return(ddate)
}

celfileChip <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	return(as.character(h$cdfName))
}

## wrapper for jetset
mappingJetset <- function(probes, platform=c("hgu133plus2", "hgu95av2", "hgu133a", "u133x3p")) {
	require(jetset)
	js <- jetset::jscores(chip="hgu133plus2", probeset=probes)
	js <- js[probes, , drop=FALSE]
	## identify the best probeset for each entrez gene id
	geneid1 <- as.character(js[ ,"EntrezID"])
	names(geneid1) <- rownames(js)
	geneid2 <- sort(unique(geneid1))
	names(geneid2) <- paste("geneid", geneid2, sep=".")
	gix1 <- !is.na(geneid1)
	gix2 <- !is.na(geneid2)
	geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
	## probes corresponding to common gene ids
	gg <- names(geneid1)[is.element(geneid1, geneid.common)]
	gid <- geneid1[is.element(geneid1, geneid.common)]
	## duplicated gene ids
	gid.dupl <- unique(gid[duplicated(gid)])
	gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
	## unique gene ids
	gid.uniq <- gid[!is.element(gid, gid.dupl)]
	gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
	## which are the best probe for each gene
	js <- data.frame(js, "best"=FALSE)
	js[gg.uniq, "best"] <- TRUE
	## data for duplicated gene ids
	if(length(gid.dupl) > 0) {
		library(jetset)
		## use jetset oevrall score to select the best probeset
		myscore <- js[gg.dupl,"overall"]
		myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
		myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
		myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
		js[myscore[ ,"probe"], "best"] <- TRUE
	}
	annot <- data.frame("probe"=rownames(js), "EntrezGene.ID"=js[ ,"EntrezID"], js)
	annot <- annot[probes, , drop=FALSE]
	return (annot)
}

CombineExpressionSetList <-
function(...) {
   args <- list(...)
   nargs <- length(args)
   if (nargs <= 1) {
     if (nargs == 1 && is.list(args[[1]])) {
       do.call("CombineExpressionSetList", args[[1]])
     } else {
       return (args[[1]])
     }
   } else if (nargs == 2) {
     return (Biobase::combine(args[[1]], args[[2]]))
   } else {
     return (Biobase::combine(args[[1]], CombineExpressionSetList(args[-1])))
   }
}

intersectList <-
function(...) {
   args <- list(...)
   nargs <- length(args)
   if (nargs <= 1) {
     if (nargs == 1 && is.list(args[[1]])) {
       do.call("intersectList", args[[1]])
     } else {
       return (args[[1]])
     }
   } else if (nargs == 2) {
     return (intersect(args[[1]], args[[2]]))
   } else {
     return (intersect(args[[1]], intersectList(args[-1])))
   }
}


## Normalize rows of a matrix to have the same quantiles, allowing for missing values.
## inpired from Gordon Smyth (normalizaeQuantiles from the limma package)
##
## input:
##  A: numeric matrix. Missing values are allowed.
##  ties: logical. If ‘TRUE’, ties in each row of ‘A’ are treated in careful way. tied values will be normalized to the mean of the corresponding pooled quantiles.
## normvector: numeric vector of values corresponding to the quantiles distribution to fit. Note that names of 'normvector' should be properly set to the colnames of 'A'
##
## output: quantil normalized matrix
`normQuant` <- 
function(A, ties=TRUE, normvector) {

  if(!missing(normvector)) {
    normvector <- sort(normvector, method="quick")
    if(length(normvector) != ncol(A)) {
      stop("length of normvector must be equal to the number of columns of A")
    }
  }
  
	n <- dim(A)
	if (is.null(n)) { return(A) }
	if (n[1] == 1) { return(A) }
	O <- S <- array( , n)
	nobs <- rep(n[2], n[1])
	i <- (0:(n[2]-1)) / (n[2]-1)
	for (j in 1:n[1]) {
		Si <- sort(A[j, ], method="quick", index.return=TRUE)
		nobsj <- length(Si$x)
		if(nobsj < n[2]) {
			nobs[j] <- nobsj
			isna <- is.na(A[j, ])
			S[j, ] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties="ordered")$y
			O[j, !isna] <- ((1:n[2])[!isna])[Si$ix]
		} else {
			S[j, ] <- Si$x
			O[j, ] <- Si$ix
		}
	}
	if(!missing(normvector)) {
    m <- normvector
  } else {
    m <- colMeans(S)
  }
	for (j in 1:n[1]) {
		if (ties) { r <- rank(A[j, ]) }
		if (nobs[j] < n[2]) {
			isna <- is.na(A[j, ])
			if (ties) {
				A[j, !isna] <- approx(i, m, (r[!isna]-1)/(nobs[j]-1), ties="ordered")$y
      } else {
				A[j, O[j, !isna]] <- approx(i, m, (0:(nobs[j]-1))/(nobs[j]-1), ties="ordered")$y
      }
		} else {
			if(ties)
				A[j, ] <- approx(i, m, (r-1)/(n[2]-1), ties="ordered")$y
			else
				A[j, O[j, ]] <- m
		}
	}
	return(A)
}

## unsupervised feature selection using sparse pca

featureSelectionSPCA <- function (exprs.obj, K=10, max.var=0.8, para=1e3) {

  spc <- elasticnet::arrayspc(x=t(exprs(exprs.obj)), K=K, para=para)
  ## choose the number of pc to either summarize as much as max.var or K pc otherwise
  
  myx <- abs(spc1$loadings) > 1e-2
  # myx <- order(abs(spc1$loadings), decreasing=TRUE)[1:500]
  genes <- featureNames(exprs.obj)[myx]
  
}
