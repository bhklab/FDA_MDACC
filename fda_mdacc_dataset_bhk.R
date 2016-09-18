########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 15, 2015
########################

########################
## normalization using frma
########################

## dataset.fn is the file name to save the normalized data
if(!file.exists(dataset.fn)) {

  library(affy)
  library(affyio)
  library(Hmisc)
  library(frma)
  
  ## CEL file names
  celfn <- list.celfiles(rawpath, full.names=TRUE)
  celfns <- list.celfiles(rawpath, full.names=FALSE)
  ## experiments' names
  names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
  ## chip type and date
  chipt <- sapply(celfn, celfileChip)
  chipd <- t(sapply(celfn, celfileDateHour))
  uchipt <- sort(unique(chipt))
  data.ge <- lapply(uchipt, function (x, y, z, w) {
    iix <- NULL
    switch (x,
      "HG-U133A" = {
        ## chip hgu133a
        library(hgu133afrmavecs)
        data(hgu133afrmavecs)
        library(hgu133acdf)
        data(hgu133acdf)
      },
      "HG-U133_Plus_2" = {
        ## chip hgu133plus2
        library(hgu133plus2frmavecs)
        data(hgu133plus2frmavecs)
        library(hgu133plus2cdf)
        data(hgu133plus2cdf)
      }
    )
    ## platform name for jetset
    platf <- tolower(gsub(badchars, "", x))
    ## CEL file names
    iix <- w[!is.na(y) & y == x]
    ## reorder CEL files by hybridization time
    iix <- iix[order(z[names(iix), "day"], z[names(iix), "hour"])]
    
    message(sprintf("Normalize gene expression data for chip %s", x))
    ## RMA
    # rr <- just.rma(filenames=iix)
    # FRMA normalization using parallel
    splitix <- parallel::splitIndices(nx=length(iix), ncl=nbcore)
    splitix <- splitix[sapply(splitix, length) > 0]
    res <- parallel::mclapply(splitix, function(x, w) {
      tt <- w[x]
      names(tt) <- NULL
      abatch <- affy::read.affybatch(filenames=tt)
      rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE)
      ## rename sample names
      colnames(exprs(rr)) <- rownames(pData(rr)) <- gsub("[.]CEL$", "", gsub("[.]gz$", "", colnames(exprs(rr))))
      return (rr)
    }, w=iix)
     res2 <- CombineExpressionSetList(res)

     ## build annotation matrix
     message("Build annotation matrix")
     annot <- mappingJetset(probes=rownames(exprs(res2)), platform=platf)
     fData(res2) <- annot
     return (res2)
  }, y=chipt, z=chipd, w=celfn)
  names(data.ge) <- uchipt
  save(list=c("data.ge"), compress=TRUE, file=dataset.fn)
} else {
  load(dataset.fn)
}



myfn <- sprintf("%s_merged.%s", file_path_sans_ext(dataset.fn), file_ext(dataset.fn))
if(!file.exists(myfn)) {
  ## subset to common probes
  gl <- lapply(data.ge, function (x) { return (featureNames(x)) })
  gl <- intersectList(gl)
  gl <- fData(data.ge[["HG-U133A"]])[gl, , drop=FALSE]
  data.ge2 <- lapply(data.ge, function (x, y) {
    exprs(x) <- exprs(x)[rownames(gl), , drop=FALSE]
    fData(x) <- gl
    annotation(x) <- "hgu133a"
    return (x)
  }, y=gl)

  ## quantile normalization between hgu133a and hgu133plus2
  tt <- apply(Biobase::exprs(data.ge2[["HG-U133A"]]), 2, sort, method="quick")
  qnormvec.hgu133a <- apply(tt, 1, median, na.rm=TRUE)
  ## normalize the data
  dd <- lapply(data.ge2, function (x, y) {
    res <- t(normQuant(A=t(exprs(x)), ties=TRUE, normvector=y))
    return (res)
  }, y=qnormvec.hgu133a)
  dd <- do.call(cbind, dd)

  data.ge3 <- CombineExpressionSetList(data.ge2)
  exprs(data.ge3) <- dd[rownames(exprs(data.ge3)), colnames(exprs(data.ge3))]
  save(list=c("data.ge3"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

########################
## clustering
########################

myfn <- sprintf("%s_merged_clustering.%s", file_path_sans_ext(dataset.fn), file_ext(dataset.fn))
if(!file.exists(myfn)) {
  hcl <- pvcl <- NULL
  
  ## select genes using sparse pca
  spc1 <- elasticnet::arrayspc(x=t(exprs(data.ge3)), K=1, para=1e4)
  myx <- abs(spc1$loadings) > 1e-2
  # myx <- order(abs(spc1$loadings), decreasing=TRUE)[1:500]
  genes <- featureNames(data.ge3)[myx]

  ## clustering using hclust
  hcl <- amap::hcluster(x=t(exprs(data.ge3)[genes, , drop=FALSE]), link="complete", method="correlation")
  hcl.k2 <- cutree(hcl, k=2)
  table("CHIPTYPE"=chipt[names(hcl.k2)], "HCLUST2"=hcl.k2)
  
  pdf(file.path(saveres, "hclust_spca_ge_all.pdf"), width=35)
  plot(hcl)
  dev.off()
  
  ## clustering using pvclust
  pvcl <- pvclust::pvclust(data=exprs(data.ge3)[genes, , drop=FALSE], method.hclust="complete", method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

  pdf(file.path(saveres, "hclust_spca_ge_all.pdf"), width=35)
  plot(pvcl)
  dev.off()

  save(list=c("hcl", "pvcl"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}


## end
  


