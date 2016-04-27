## install the new release version of PharmacoGx (version 1.1.4)
if(!(require(PharmacoGx))){
  source('http://bioconductor.org/biocLite.R')
  biocLite('PharmacoGx')
} 

library(PharmacoGx)
library(Biobase)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 4
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

## directory where all the analysis results will be stored
saveres <- "Output"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }

#################################################
### functions

myScatterPlot <- function(x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...) {
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  if (length(col) != length(x)) {
    col <- rep(col, length.out=length(x))
  }
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  col <- col[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             rr <- plot(x=x, y=y, col=col, pch=pch, ...)
           },
           "transparent"={
             myrgb <- sapply(col, grDevices::col2rgb, alpha=FALSE) / 255
             myrgb <- apply(myrgb, 2, function (x, transparency) {
               return (rgb(red=x[1], green=x[2], blue=x[3], alpha=transparency, maxColorValue=1))
             }, transparency=transparency)
             rr <- plot(x=x, y=y, pch=pch, col=myrgb, ...)
           },
           "smooth"={
             rr <- smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
           }
    )
  }
  
  invisible(rr)
}

#################################################
###  Nature 2013 PSets, the old version of datasets which have been used in Haibe-Kains et. al.
gdsc.2013 <- downloadPSet("GDSC_2013")
ccle.2013 <- downloadPSet("CCLE_2013")

### read old curation file which contains the curated cellid for 471 cell lines
##Data_File_2.csv contains the old curation has been done by Haibe-Kains et. al. and has been used in both inconsistency/agreement papers
common.cells.2013 <- read.csv("Data_File_2.csv", stringsAsFactors=FALSE)

### subset psets to old cutation(471 cell lines)
gdsc.2013 <- subsetTo(gdsc.2013, cells=common.cells.2013$cellid)
ccle.2013 <- subsetTo(ccle.2013, cells=common.cells.2013$cellid)

### find the intersection of psets
common.2013 <- intersectPSet(pSets=list("CCLE"=ccle.2013, "GDSC"=gdsc.2013),intersectOn=c("cell.lines","drugs"), strictIntersect=TRUE)
rr2 <- sensNumber(common.2013$CCLE)
rr1 <- sensNumber(common.2013$GDSC)
table(rr1 == rr2)
sensitivity.ids.2013 <- paste(sensitivityInfo(common.2013$CCLE)[, "cellid"], sensitivityInfo(common.2013$CCLE)[, "drugid"], sep="_")

### load the newest version of psets containing the most recent released data
GDSC <- downloadPSet(name="GDSC")
CCLE <- downloadPSet(name="CCLE")
### find the intersection of psets
common <- intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC),intersectOn=c("cell.lines","drugs"), strictIntersect=TRUE)

rr2 <- sensNumber(common$CCLE)
rr1 <- sensNumber(common$GDSC)
table(rr1 == rr2)

### The number of experiments missed because of using an old version of data and old curation
sensitivity.ids <- paste(sensitivityInfo(common$CCLE)[, "cellid"], sensitivityInfo(common$CCLE)[, "drugid"], sep="_")
length(setdiff(sensitivity.ids, sensitivity.ids.2013))

### outdated data
### GDSC
sensitivity.ids.2013 <- paste(sensitivityInfo(gdsc.2013)[, "cellid"], sensitivityInfo(gdsc.2013)[, "drugid"], sep="_")
sensitivity.ids <- paste(sensitivityInfo(GDSC)[, "cellid"], sensitivityInfo(GDSC)[, "drugid"], sep="_")
sensitivity.ids.common <- intersect(sensitivity.ids.2013, sensitivity.ids)
cells <- unlist(lapply(strsplit(sensitivity.ids.common, split="_"), function(x){x[[1]]}))
drugs <- unlist(lapply(strsplit(sensitivity.ids.common, split="_"), function(x){x[[2]]}))

sensitivity.ids.all <- NULL
sensitivity.ids.2013 <- NULL
for(i in 1:length(cells)){
  x <- rownames(sensitivityInfo(GDSC))[which(sensitivityInfo(GDSC)[, "cellid"] == cells[i] & sensitivityInfo(GDSC)[, "drugid"] == drugs[i])]
  y <- rownames(sensitivityInfo(gdsc.2013))[which(sensitivityInfo(gdsc.2013)[, "cellid"] == cells[i] & sensitivityInfo(gdsc.2013)[, "drugid"] == drugs[i])]
  sensitivity.ids.all <- c(sensitivity.ids.all, x)
  sensitivity.ids.2013 <- c(sensitivity.ids.2013, y)
}
xx <- GDSC@sensitivity$profiles[sensitivity.ids.all, "ic50_published"] 
yy <- gdsc.2013@sensitivity$profiles[sensitivity.ids.2013, "ic50_published"]
table(xx == yy)["FALSE"]

xx <- GDSC@sensitivity$profiles[sensitivity.ids.all, "auc_published"] 
yy <- gdsc.2013@sensitivity$profiles[sensitivity.ids.2013, "auc_published"]
table(xx == yy)

### CCLE
sensitivity.ids.2013 <- paste(sensitivityInfo(ccle.2013)[, "cellid"], sensitivityInfo(ccle.2013)[, "drugid"], sep="_")
sensitivity.ids <- paste(sensitivityInfo(CCLE)[, "cellid"], sensitivityInfo(CCLE)[, "drugid"], sep="_")
sensitivity.ids.common <- intersect(sensitivity.ids.2013, sensitivity.ids)
cells <- unlist(lapply(strsplit(sensitivity.ids.common, split="_"), function(x){x[[1]]}))
drugs <- unlist(lapply(strsplit(sensitivity.ids.common, split="_"), function(x){x[[2]]}))

sensitivity.ids.all <- NULL
sensitivity.ids.2013 <- NULL
for(i in 1:length(cells)){
  x <- rownames(sensitivityInfo(CCLE))[which(sensitivityInfo(CCLE)[, "cellid"] == cells[i] & sensitivityInfo(CCLE)[, "drugid"] == drugs[i])]
  y <- rownames(sensitivityInfo(ccle.2013))[which(sensitivityInfo(ccle.2013)[, "cellid"] == cells[i] & sensitivityInfo(ccle.2013)[, "drugid"] == drugs[i])]
  sensitivity.ids.all <- c(sensitivity.ids.all, x)
  sensitivity.ids.2013 <- c(sensitivity.ids.2013, y)
}
xx <- CCLE@sensitivity$profiles[sensitivity.ids.all, "ic50_published"] 
yy <- ccle.2013@sensitivity$profiles[sensitivity.ids.2013, "ic50_published"]
table(xx == yy)["FALSE"]

xx <- CCLE@sensitivity$profiles[sensitivity.ids.all, "auc_published"] 
yy <- ccle.2013@sensitivity$profiles[sensitivity.ids.2013, "auc_published"]
table(xx == yy)

###  consistency of tehcnical replicates in GDSC sites
### AZD6482
ids <- unlist(strsplit(GDSC@drug[which(rownames(GDSC@drug) == "AZD6482"), "drugid"], split="/"))

tt <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "ic50_published"]
ss <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)),"ic50_published"]

names(tt) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "cellid"]
names(ss) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)), "cellid"]

nn <- intersect(names(tt), names(ss))
tt <- -log10(tt)
ss <- -log10(ss)

length(nn)
cc <- cor.test(tt[nn], ss[nn])
ma <- max(tt[nn], ss[nn])
mi <- min(tt[nn], ss[nn])
ll <- lm(ss[nn] ~ tt[nn])

##Extended Data Figure 1.b: consistency of sensitivity prfiles across replicated experiments in GDSC: AZD6482
pdf(file.path(saveres, "AZD6482_sites_replicates.pdf"), height=5, width=5)
myScatterPlot(tt[nn], ss[nn], xlab="-Log10 IC50(WTSI)", ylab="-Log10 IC50(MGH)", xlim=c(mi, ma), ylim=c(mi, ma))
abline(0, 1, col="black", lty=2)
abline(coef(ll)[1], coef(ll)[2], col="gray")

legend("topleft", legend=sprintf("PCC=%s, p=%s\n#cell lines=%s", round(cc$estimate, digits=2), cc$p.value, length(nn)), bty="n")
dev.off()
#####################################

### Camptothecin
message("Download drug information")
if(!file.exists(file.path(saveres, "nature_supplementary_information.xls"))){
  dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(saveres, "nature11005-s2.zip"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  ff <- as.character(unzip(zipfile=file.path(saveres, "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(saveres, "nature11005-s2.zip"), exdir=saveres)
  file.copy(from=file.path(saveres, ff), to=file.path(saveres, "nature_supplementary_information.xls"))
}

if(!file.exists(file.path(saveres, "drugpheno.nature.RData"))) {
  drugpheno.nature <- gdata::read.xls(xls=file.path(saveres, "nature_supplementary_information.xls"), sheet=2)
  drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
  save(drugpheno.nature, file=file.path(saveres, "drugpheno.nature.RData"))
} else {
  load(file.path(saveres, "drugpheno.nature.RData"))
}

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
coln2 <- gsub(" ", "", sapply(drugpheno.nature[1,], as.character))
coln2[coln2 == ""] <- NA
drugpheno.nature <- drugpheno.nature[-1, ,drop=FALSE]
coln <- colnames(drugpheno.nature)
coln2[is.na(coln2)] <- coln[is.na(coln2)]
coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
colnames(drugpheno.nature) <- coln2
myx <- sapply(strsplit(colnames(drugpheno.nature), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
ic50 <- drugpheno.nature[ , myx, drop=FALSE]
nn <- dimnames(ic50)
nn[[2]] <- gsub("_IC_50", "", nn[[2]])
ic50 <- apply(ic50, 2, as.numeric)
dimnames(ic50) <- nn
ic50 <- exp(ic50) / 10^6
###  camptothecin
yy <- -log10(ic50[ , "drugid_195", drop=FALSE])
xx <- -log10(ic50[ , "drugid_1003", drop=FALSE])
ccix <- complete.cases(xx, yy)
nnn <- sum(ccix)
cc <- cor.test(x=xx, y=yy, method="pearson", use="complete.obs")
ma <- max(xx, yy, na.rm=T)
mi <- min(xx, yy, na.rm=T)
ll <- lm(yy ~ xx)

##Extended Data Figure 1.a: consistency of sensitivity prfiles across replicated experiments in GDSC: Camptothecin
pdf(file.path(saveres, "Camptothecin_sites_replicates.pdf"), height=5, width=5)
myScatterPlot(xx, yy, xlab="-Log10 IC50(WTSI)", ylab="-Log10 IC50(MGH)", xlim=c(mi, ma), ylim=c(mi, ma))
abline(0, 1, col="black", lty=2)
abline(coef(ll)[1],  coef(ll)[2], col="gray")
legend("topleft", legend=sprintf("PCC=%s, p=%s\n#cell lines=%s", round(cc$estimate, digits=2), cc$p.value, nnn), bty="n")
dev.off()
#####################################

###  consistency of molecular profiles versus sensitivity profiles

### consistency of gene expression
rna.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common$GDSC, mDataType="rna", summary.stat="median"))
rna.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common$CCLE, mDataType="rna", summary.stat="median"))
###  common genes
gix <- intersect(rownames(rna.gdsc.common), rownames(rna.ccle.common))
rna.gdsc.common <- rna.gdsc.common[gix, ]
rna.ccle.common <- rna.ccle.common[gix, ]

myfn <- file.path(saveres, sprintf("rna_consistency.RData"))
if (!file.exists(myfn)) {
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  consis <- matrix(NA, nrow=length(gix), ncol=2, dimnames=list(gix,c("estimate", "p")))
  ii <- 1
  for(i in gix){
    res <- try(cor.test(rna.gdsc.common[i, ], rna.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, ] <- c(res$estimate, res$p.value)
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myfn)
}else {
  load(myfn)
}
rna.consistency <- consis

### consistency of cnv profiles
cnv.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common$GDSC, mDataType="cnv", summary.stat="median"))
cnv.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common$CCLE, mDataType="cnv", summary.stat="median"))
###  common genes
gix <- intersect(rownames(cnv.gdsc.common), rownames(cnv.ccle.common))
cnv.gdsc.common <- cnv.gdsc.common[gix, ]
cnv.ccle.common <- cnv.ccle.common[gix, ]

myfn <- file.path(saveres, sprintf("cnv_consistency.RData"))
if (!file.exists(myfn)) {
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  consis <- matrix(NA, nrow=length(gix), ncol=2, dimnames=list(gix,c("estimate", "p")))
  ii <- 1
  for(i in gix){
    res <- try(cor.test(cnv.gdsc.common[i, ], cnv.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, ] <- c(res$estimate, res$p.value)
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myfn)
}else {
  load(myfn)
}
cnv.consistency <- consis

### consistency of mutation profiles
mutation.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common$GDSC, mDataType="mutation", summary.stat="or"))
mutation.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common$CCLE, mDataType="mutation", summary.stat="or"))
### common genes
gix <- intersect(rownames(mutation.gdsc.common), rownames(mutation.ccle.common))
mutation.gdsc.common <- mutation.gdsc.common[gix, ]
mutation.ccle.common <- mutation.ccle.common[gix, ]

consis <- matrix(NA, nrow=length(gix), ncol=2, dimnames=list(gix,c("estimate", "p")))
for(i in gix){
  mutation.gdsc.common.bin <- factor(ifelse(mutation.gdsc.common[i, ] == 1, "mutated", "wild.type"), levels=c("mutated","wild.type"))
  mutation.ccle.common.bin <- factor(ifelse(mutation.ccle.common[i, ] == 1, "mutated", "wild.type"), levels=c("mutated","wild.type"))
  
  ###  kappa
  tt <- table(x=mutation.gdsc.common.bin, y=mutation.ccle.common.bin)
  res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
  if (class(res) != "try-error") {
    consis[i, "estimate"] <- res$kappa
  }
}
mutation.consistency <- consis



### consistency of IC50 profiles
ic50.gdsc.common <- summarizeSensitivityProfiles(pSet=common$GDSC, sensitivity.measure="ic50_published", summary.stat="median")
ic50.ccle.common <- summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="ic50_published", summary.stat="median")
drugs <- drugNames(common$GDSC)

max.conc.ccle.common <- max.conc.gdsc.common <- matrix(NA, nrow=nrow(ic50.gdsc.common), ncol=ncol(ic50.gdsc.common), dimnames=dimnames(ic50.gdsc.common))
for(cellid in colnames(max.conc.ccle.common)) {
  for(drugid in rownames(max.conc.ccle.common)) {
    tt <- max(as.numeric(sensitivityInfo(common$CCLE)[which(sensitivityInfo(common$CCLE)[,"cellid"] == cellid & sensitivityInfo(common$CCLE)[,"drugid"] == drugid), grep("Dose", colnames(sensitivityInfo(common$CCLE)))]))
    max.conc.ccle.common[drugid, cellid] <- ifelse(length(tt) > 0, tt, NA)
    tt <- as.numeric(sensitivityInfo(common$GDSC)[which(sensitivityInfo(common$GDSC)[,"cellid"] == cellid & sensitivityInfo(common$GDSC)[,"drugid"] == drugid), "max.Dose.uM"])
    max.conc.gdsc.common[drugid, cellid] <- ifelse(length(tt) > 0, tt, NA)
  }
}

### Cutoff set to 1uM
cc <- 1
consis <- array(NA, dim=c(length(drugs), 2, 2), dimnames=list(drugs, c("PCC", "kappa"), c("estimate", "p")))
for(i in drugs){
  ic50.gdsc.common.i <- ic50.gdsc.common[i, ]
  max.gdsc.common.i <- max.conc.gdsc.common[i, ]
  tt <- ic50.gdsc.common.i > max.gdsc.common.i
  ic50.gdsc.common.i[which(tt == T)] <- max.gdsc.common.i[which(tt == T)]
  
  ic50.ccle.common.i <- ic50.ccle.common[i, ]
  max.ccle.common.i <- max.conc.ccle.common[i, ]
  tt <- ic50.ccle.common.i > max.ccle.common.i
  ic50.ccle.common.i[which(tt == T)] <- max.ccle.common.i[which(tt == T)]
  
  ### PCC
  res <- try(cor.test(ic50.gdsc.common.i, ic50.ccle.common.i, method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
  if (class(res) != "try-error") {
    consis[i, "PCC", ] <- c(res$estimate, res$p.value)
  }
  
  ic50.gdsc.common.bin <- factor(ifelse(ic50.gdsc.common[i, ] < cc, "sensitive", "resistent"), levels=c("sensitive", "resistent"))
  ic50.ccle.common.bin <- factor(ifelse(ic50.ccle.common[i, ] < cc, "sensitive", "resistent"), levels=c("sensitive", "resistent"))
  
  ###  kappa
  tt <- table(x=ic50.gdsc.common.bin, y=ic50.ccle.common.bin)
  res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
  if (class(res) != "try-error") {
    consis[i, "kappa", "estimate"] <- res$kappa
  }
  
}
ic50.consistency <- consis

### AUC cutoff
gdsc.auc <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=common$GDSC, sensitivity.measure="auc_published", summary="median"))
ccle.auc <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="auc_published", summary="median"))

cols <- matrix("#D7191C", ncol=ncol(auc.gdsc.common), nrow=nrow(auc.gdsc.common), dimnames=dimnames(auc.gdsc.common))
cols[which(ic50.ccle.common < cc & ic50.gdsc.common < cc)] <- "#2B83BA"
cols[which(ic50.ccle.common < cc & ic50.gdsc.common >= cc)] <- "#FDAE61"
cols[which(ic50.ccle.common >= cc & ic50.gdsc.common < cc)] <- "#ABDDA4"

pdf(file.path(saveres, "auc_colored_ic50_cutoff_1uM.pdf"), height=5, width=5)
myScatterPlot(x=gdsc.auc, y=ccle.auc, method="transparent", transparency=0.70, smooth.pch=".", pch=16, minp=50, col=t(cols), xlab="AUC (GDSC)", ylab="AUC (CCLE)", ylim=c(0, 1), xlim=c(0, 1))
dev.off()
### consistency of AUC profiles
cc <- 0.2
auc.gdsc.common <- summarizeSensitivityProfiles(pSet=common$GDSC, sensitivity.measure="auc_published", summary.stat="median")
auc.ccle.common <- summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="auc_published", summary.stat="median")
drugs <- drugNames(common$GDSC)

consis <- array(NA, dim=c(length(drugs), 2, 2), dimnames=list(drugs, c("PCC", "kappa"), c("estimate", "p")))
for(i in drugs){
  ### PCC
  res <- try(cor.test(auc.gdsc.common[i, ], auc.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
  if (class(res) != "try-error") {
    consis[i, "PCC", ] <- c(res$estimate, res$p.value)
  }
  
  auc.gdsc.common.bin <- factor(ifelse(auc.gdsc.common[i, ] < cc, "sensitive", "resistent"), levels=c("sensitive", "resistent"))
  auc.ccle.common.bin <- factor(ifelse(auc.ccle.common[i, ] < cc, "sensitive", "resistent"), levels=c("sensitive", "resistent"))
  
  ###  kappa
  tt <- table(x=auc.gdsc.common.bin, y=auc.ccle.common.bin)
  res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
  if (class(res) != "try-error") {
    consis[i, "kappa", "estimate"] <- res$kappa
  }
  
}
auc.consistency <- consis

### Extended Data Figure 1.c and 1.d: consistency of molecular and pharmacological profiles using continuous and binarized values for cnv, expression, auc and ic50
pdf(file.path(saveres, "boxplot_molecular_sensitivity_consistency.pdf"), height=5, width=10)
par(mfrow=c(1,2))

### PCC
xx <- list(
  "cnv"=cnv.consistency[ , "estimate"],
  "expression"=rna.consistency[ ,"estimate"],
  "auc"=auc.consistency[ , "PCC", "estimate"],
  "ic50"=ic50.consistency[ , "PCC", "estimate"]
)
yylim <- c(-0.1, 1)
bp <- boxplot(xx, ylim=yylim, col="lightgrey", border="black", outline=FALSE, pars=list(outcol="black", outpch=20, outcex=0.5), ylab="PCC", xaxt="n")
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=1:length(xx), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2, cex=0.8)
### kappa
xx <- list(
  "mutation"=mutation.consistency[ ,"estimate"],
  "auc"=auc.consistency[ , "kappa", "estimate"],
  "ic50"=ic50.consistency[ , "kappa", "estimate"]
)
yylim <- c(-0.1, 1)
bp <- boxplot(xx, ylim=yylim, col="lightgrey", border="black", outline=FALSE, pars=list(outcol="black", outpch=20, outcex=0.5), ylab="kappa", xaxt="n")
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=1:length(xx), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2, cex=0.8)
dev.off()
###########################################

writeLines(utils::toLatex(sessionInfo()), file.path(saveres, "sessionInfo.tex"))


## end

