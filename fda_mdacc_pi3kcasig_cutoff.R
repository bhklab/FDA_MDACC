########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 11, 2017
########################

load("loi2008.RData")
ge <- data
ge.annot <- annot
clin <- demo

### subtyping
sbt <- molecular.subtyping(sbt.model="scmod1", data=ge, annot=ge.annot, do.mapping=FALSE)
iix <- !is.na(sbt$subtype) & (is.element(sbt$subtype, c("ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")))

### select the ER+/HER2- patients
iix.wt <- rownames(clin)[complete.cases(clin[ , c("PI3KCA.mutation.exon20.status")]) & clin[ , "PI3KCA.mutation.exon20.status"] == 0]
iix.mut <- rownames(clin)[complete.cases(clin[ , c("PI3KCA.mutation.exon20.status")]) & clin[ , "PI3KCA.mutation.exon20.status"] == 1]

sig$score <- genefu::pik3cags(data=ge[iix, , drop=FALSE], annot=ge.annot, do.mapping=FALSE)
### rescale
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
### define cutoff
m.wt <- mean(sig$score[iix.wt], na.rm = TRUE)
m.mut <- mean(sig$score[iix.mut], na.rm = TRUE)
pi3kcasig.lum.rescale.cutoff <- m.wt + (m.mut - m.wt)/2

