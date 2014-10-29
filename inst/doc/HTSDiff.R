### R code from vignette source 'HTSDiff.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(digits=3, width=100)


###################################################
### code chunk number 2: quickstart (eval = FALSE)
###################################################
## y <- read.table("counts.txt")
## conds <- c(1,1,2,2)
## mod <- HTSDiff(y, conds)
## DEresults <- mod$res


###################################################
### code chunk number 3: loaddat
###################################################
library(HTSDiff)
data(initialDataset)
y <- initialDataset[,c("BF1", "BF2", "F1", "F2")]
## Fix gene IDs as row names
rownames(y) <- initialDataset[,1]
head(y)
conds <- c("BF","BF","F","F")
y <- y[rowSums(y)>0,]
dim(y)


###################################################
### code chunk number 4: diffanalysis
###################################################
set.seed(12345)
DEtest <- HTSDiff(counts=y, conds=conds)


###################################################
### code chunk number 5: res
###################################################
res <- DEtest$res
head(res)


###################################################
### code chunk number 6: DE
###################################################
table(res$DE)


###################################################
### code chunk number 7: PMM
###################################################
names(DEtest$PMM)


###################################################
### code chunk number 8: summary
###################################################
summary(DEtest$PMM)


###################################################
### code chunk number 9: edgeR
###################################################
library(edgeR)
y <- DGEList(counts=y, group=conds)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
de <- decideTestsDGE(et, p=0.05, adjust="BH")
summary(de)


###################################################
### code chunk number 10: venn
###################################################
tab <- table(abs(de),res$DE)
rownames(tab) <- c("NDE (edgeR)", "DE (edgeR)")
colnames(tab) <- c("NDE (HTSDiff)", "DE (HTSDiff)")
tab


###################################################
### code chunk number 11: bcv
###################################################

A <- y$AveLogCPM
disp <- getDispersion(y)
colors <- ifelse(abs(de)==1 & res$DE==TRUE, "grey5", "grey70")
colors <- ifelse(abs(de)==1 & res$DE==FALSE, "blue", colors)
colors <- ifelse(abs(de)==0 & res$DE==TRUE, "red", colors)
cex <- rep(0.4, length(abs(de)))
cex <- ifelse(abs(de)==1 & res$DE==FALSE, 0.7, cex)
cex <- ifelse(abs(de)==0 & res$DE==TRUE, 0.7, cex)
plot(A, sqrt(disp), xlab="Average log CPM", ylab="Biological coef of variation", 
     type="n")
points(A, sqrt(y$tagwise.dispersion), pch = 16, cex=cex, col=colors)
legend("topright", c("NDE", "DE (HTSDiff & edgeR)", "DE (edgeR)", "DE (HTSDiff"), 
       col=c("grey70","grey5","blue", "red"), bty="n", pch=16, cex=0.7)



###################################################
### code chunk number 12: sessi
###################################################
sessionInfo()


