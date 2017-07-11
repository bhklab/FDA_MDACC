########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 11, 2017
########################

load("loi2008.RData")
ge <- data
ge.annot <- annot
clin <- demo

### ER status by mRNA
er.status <- genefu::bimod(x=mod1$ESR1[1, ], data=ge, annot=ge.annot, do.mapping=FALSE, model="V", do.scale=FALSE)$status

### select the ER+ grade 1 patients
iix1 <- rownames(clin)[complete.cases(clin[ , c("er", "grade")]) & er.status == 1 & clin[ , "grade"] == 1]
iix3 <- rownames(clin)[complete.cases(clin[ , c("er", "grade")]) & er.status == 1 & clin[ , "grade"] == 3]

sig <- ggi(data=ge[er.status == 1, , drop=FALSE], annot=ge.annot, mapping=FALSE)
### rescale
sig$score <- genefu::rescale(sig$score, na.rm=TRUE, q=0.05)
### define cutoff
mhg1 <- mean(sig$score[iix1], na.rm = TRUE)
mhg3 <- mean(sig$score[iix3], na.rm = TRUE)
ggi.erp.rescale.cutoff <- mhg1 + (mhg3 - mhg1)/2


