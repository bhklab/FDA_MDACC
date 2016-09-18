########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 15, 2015
########################

## data.ge3 contains the normalized data

colo <- c("darkblue", "darkorange", "darkred")

res.scores <- res.risks <- data.frame(matrix(paste(sampleNames(data.ge3), "CEL", sep="."), nrow=length(sampleNames(data.ge3)), ncol=1, dimnames=list(sampleNames(data.ge3), "CEL_File_Name")))

#################################################
## signatures based on affymetrix data

##############
## GGI
##############
sig.name <- "GGI"
## without histological grade, only continuous score can be computed
sig <- genefu::ggi(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=FALSE)
## rescale value between 0 and 1
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
sig$score[!is.na(sig$score) & sig$score < 0] <- 0
sig$score[!is.na(sig$score) & sig$score > 1] <- 1

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## GENIUS
##############
sig.name <- "GENIUS"
sig <- genefu::genius(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=FALSE)
## rescale value between 0 and 1
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
sig$score[!is.na(sig$score) & sig$score < 0] <- 0
sig$score[!is.na(sig$score) & sig$score > 1] <- 1

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## TAMR13
##############
sig.name <- "TAMR13"
sig <- genefu::tamr13(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=FALSE)
## rescale value between 0 and 1
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
sig$score[!is.na(sig$score) & sig$score < 0] <- 0
sig$score[!is.na(sig$score) & sig$score > 1] <- 1

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## PIK3CASIG
##############
sig.name <- "PIK3CASIG"
sig$score <- genefu::pik3cags(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=FALSE)
sig$risk <- rep(NA, length(sig$score))
names(sig$risk) <- names(sig$score)
## rescale value between 0 and 1
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
sig$score[!is.na(sig$score) & sig$score < 0] <- 0
sig$score[!is.na(sig$score) & sig$score > 1] <- 1

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL


########################
## Gene modules from Desmedt et al, CCR 2008
########################

for (i in 1:length(genefu::mod1)) {
  sig.name <- paste(names(genefu::mod1)[i], "MODULE", sep="_")
  sig$score <- genefu::sig.score(x=genefu::mod1[[i]], data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=FALSE)$score
  sig$risk <- rep(NA, length(sig$score))
  names(sig$risk) <- names(sig$score)
  ## rescale value between 0 and 1
  sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
  sig$score[!is.na(sig$score) & sig$score < 0] <- 0
  sig$score[!is.na(sig$score) & sig$score > 1] <- 1

  res.scores <- cbind(res.scores, NA)
  res.scores[names(sig$score), ncol(res.scores)] <- sig$score
  res.risks <- cbind(res.risks, NA)
  res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
  colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
  sig <- NULL
}

#################################################
## signatures developed on other platforms

########################
## Gene modules in Ignatiadis et al, JCO 2012
########################

myfn <- file.path(saveres, "ignatiadis2002_gene_modules.RData")
if(!file.exists(myfn)) {
  dwl.status <- download.file(url="http://jco.ascopubs.org/content/suppl/2012/04/16/JCO.2011.39.5624.DC1/Appendix_module_composition_file.xls", destfile=file.path(saveres, "gene_modules.xls"))
  if(dwl.status != 0) { stop("Download failed, please rerun the script!") }
  ## read xls file
  mfile <- gdata::read.xls(file.path(saveres, "gene_modules.xls"), skip=1)
  ix.delim <- c(which(mfile[ ,1] != "")[-1]-1, nrow(mfile) + 1)
  ix.f <- ix.l <- 1
  groups <- NULL
  npp <- np <- NULL
  for (i in 1:length(ix.delim)) {
  	ix.l <- ix.delim[i]
  	np <- c(np, as.character(mfile[ix.f,1]))
  	groups <- c(groups, rep(i, ix.l - ix.f + 1))
  	npp <- rbind(npp, mfile[ix.f:ix.l,2:ncol(mfile)])
  	ix.f <- ix.l + 1
  }
  ugroups <- unique(groups)
  obj <- NULL
  for (j in 1:length(ugroups)) {
  	obj <- c(obj, list(npp[groups == ugroups[j], ]))
  }
  names(obj) <- np
  ## clean modules
  gene.modules <- lapply(obj, function(x) {
    x <- x[complete.cases(x), , drop=FALSE]
    x <- x[!duplicated(x[ , "EntrezGene.ID"]), , drop=FALSE]
    rownames(x) <- paste("geneid", as.character(x[ , "EntrezGene.ID"]), sep=".")
    x <- cbind("probe"=rownames(x), x)
    return(x)
  })
  ## select "CIN70", "IRMODULE", "RAS", "MAPK", "PTEN", "AKTmTOR", "IGF1", "SRC", "MYC", "E2F3", and "BetaCatenin"
  gene.modules <- gene.modules[c("CIN70", "Immune1", "RAS", "MAPK", "PTEN", "AKTmTOR", "IGF1", "SRC", "MYC", "E2F3", "BetaCatenin")]
  names(gene.modules)[names(gene.modules) %in% "Immune1"] <- "IRMODULE"
  save(list=c("gene.modules"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

for (i in 1:length(gene.modules)) {
  sig.name <- paste(names(gene.modules)[i], "MODULE", sep="_")
  sig$score <- genefu::sig.score(x=gene.modules[[i]], data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=TRUE)$score
  sig$risk <- rep(NA, length(sig$score))
  names(sig$risk) <- names(sig$score)
  ## rescale value between 0 and 1
  sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
  sig$score[!is.na(sig$score) & sig$score < 0] <- 0
  sig$score[!is.na(sig$score) & sig$score > 1] <- 1
  
  res.scores <- cbind(res.scores, NA)
  res.scores[names(sig$score), ncol(res.scores)] <- sig$score
  res.risks <- cbind(res.risks, NA)
  res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
  colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
  sig <- NULL
}

##############
## GENE70
##############
sig.name <- "GENE70"
sig <- genefu::gene70(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=TRUE)
## rescale value between 0 and 1
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
sig$score[!is.na(sig$score) & sig$score < 0] <- 0
sig$score[!is.na(sig$score) & sig$score > 1] <- 1

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## ONCOTYPEDX
##############
sig.name <- "ONCOTYPEDX"
sig.oncotypedx.save <- genefu::sig.oncotypedx
sig.oncotypedx["GSTM1", "EntrezGene.ID"] <- "2948"
sig <- genefu::oncotypedx(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=TRUE)
## rescale value between 0 and 1
sig$score <- sig$score / 100

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## RORS
##############
sig.name <- "RORS"
sig <- genefu::rorS(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=TRUE)
## rescale value between 0 and 1
sig$score <- sig$score / 100

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL

##############
## EndoPredict
##############
sig.name <- "ENDOPREDICT"
sig <- genefu::endoPredict(data=t(exprs(data.ge3)), annot=fData(data.ge3), do.mapping=TRUE)
## rescale value between 0 and 1
sig$score <- sig$score / 15

res.scores <- cbind(res.scores, NA)
res.scores[names(sig$score), ncol(res.scores)] <- sig$score
res.risks <- cbind(res.risks, NA)
res.risks[names(sig$risk), ncol(res.risks)] <- sig$risk
colnames(res.scores)[ncol(res.scores)] <- colnames(res.risks)[ncol(res.risks)] <- sig.name
sig <- NULL


#################################################
## save results

res.all <- NULL
for (i in 2:length(res.scores)) {
  rr <- data.frame(res.scores[ , 1], res.scores[ , i], res.risks[ , i])
  dimnames(rr) <- list(rownames(res.scores), c("CEL_File_Name", "Prediction_Binary", "Prediction_Continuous"))
  res.all <- c(res.all, list(rr))
}
names(res.all) <- colnames(res.scores)[-1]
## shorten signatures names
# names(res.all) <- gsub("_MODULE$", "_M", names(res.all))

WriteXLS::WriteXLS("res.all", ExcelFileName=file.path(saveres, "fda_mdacc_signatures.xlsx"))

## end
