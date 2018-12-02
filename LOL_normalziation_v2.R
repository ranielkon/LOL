
#

#
# Run the script from the dir in which the input files are located:
#    TNF_rpkm_logR_GSE64233.txt 
#    Hs_Gencode.v25_Tx.AP.Lens.txt
#
# TNF.data.LOL.norm <- LOL.norm(exprs.file.logR="TNF_rpkm_logR_GSE64233.txt", Ylim=c(-2.5,2.5), Xlim=c(8,15))
#

#-----------------------------------------------------
#
# apply lowess norm of exprs.mat against Tx.Len (log2). 
#
#     exprs.file.logR - File with relative expression levels (log2); 
#                      1st column contains ENSEMBL gene ids (corrspoding to IDs in Hs_Gencode.v25_Tx.AP.Lens.txt)
#                      2nd column contains gene symbol
#                      All other columns contain the expression data (to be normalized)
#
#     f - vector of the variable used for normalization. (that is, Transcript length in LOL)
#     use: names(f) = rownames(M) to set gene names to f (and allow its intersection with exprs.mat)
#
# returns normalized sub.expr.mat (containg only genes for which Tx.Length.data is available)
#         LOL-norm plots are save in PDF files (file per condition)
#

LOL.norm <- function (exprs.file.logR="TNF_rpkm_logR_GSE64233.txt", Ylim=c(-2.5,2.5), Xlim=c(8,15)) {

	print("Reading Tx.Len file ...")
	Hs.Tx.LenFile <- "Hs_Gencode.v25_Tx.AP.Lens.txt"; 
	Hs.genes.len <- read.table(file=Hs.Tx.LenFile, sep="\t", header=1, row.names=1, stringsAsFactors=F); 

	print("Reading Expression file ...")
	exprs.mat.logR <- read.table(file=exprs.file.logR, sep="\t", header=1, row.names=1, stringsAsFactors=F); 

	cGenes <- intersect(rownames(exprs.mat.logR), rownames(Hs.genes.len))
	expr.submat <- exprs.mat.logR[cGenes, ]
	Tx.lens     <- log(Hs.genes.len[cGenes, "Tx.len"],2);

		# Initialize the matrix that will contain the normalized logR levels
	norm.expr.submat <- array(dim=dim(expr.submat))
	rownames(norm.expr.submat) <- rownames(expr.submat)
	colnames(norm.expr.submat) <- colnames(expr.submat)
	norm.expr.submat[, 1] <- expr.submat[, 1]; # Sym column

	conds <- colnames(expr.submat)

	x <- array(dim=c(dim(expr.submat)[1], 2))
	x[, 1] <- Tx.lens

	 # loop over data columns; run lowess.norm on each

	for (j in 2:length(conds)) { # 1st col is SYM, hence starting from j=2
		str <- sprintf("normalizing cond %s", conds[j])
		print (str)

		x[, 2] <- expr.submat[, j]
		x.norm <- lowess.norm(x, x.lab="Tx.len", Xlim, Ylim, conds[j])
		norm.expr.submat [, j] <- x.norm[, 2]
	}

	norm.expr.submat 

}

#-----------------------------------------------------

#
# x: 2D vector - 1st col is the Tx.Len and the 2nd col is log2FC
#
# returns: 
# x.norm - 2D vector with LOL normalized logR expression levels
# 

lowess.norm <- function (x, x.lab, Xlim, Ylim, cond) {

	x.norm <- array(dim=dim(x))

          	# Calculate the smoothing lowess curve for normalziation
	Len <- as.numeric(x[, 1]);
	M <- as.numeric(x[, 2]);
	LWS.curve <- loess(M ~ Len)

            # normalize the data in x according to LWS.norm curve

	k <-  predict(LWS.curve, Len); 
	M.norm <- M - k; 

	x.norm[, 1] <- Len
	x.norm[, 2] <- M.norm

		#str <- sprintf("i=%d: y=%f --> y.norm=%f", i, x[i,2], x.norm[i,2])
		#print (str)

		# Save graphs in PDF file
	print("Saving plots in PDF file ...")
	pdf.file <- sprintf("%s_LOL_norm.pdf", cond)
	pdf(pdf.file)
	par(mfrow=c(2, 1))

	str <- sprintf("%s - before normalization", cond)
	plot(x, main=str, xlim=Xlim, ylim=Ylim, xlab=x.lab, ylab="FC (log2)", pch=20, cex=0.7, col="grey")
	lines(c(0.8*Xlim[1], 1.2*Xlim[2]), c(0,0), col=4, lwd=1)
	lines(LWS.curve$x[order(LWS.curve$x)], LWS.curve$fitted[order(LWS.curve$x)], col=2, lwd=2)

	LWS.curve <- loess(as.numeric(x.norm[, 2]) ~ as.numeric(x.norm[, 1]) )
	str <- sprintf("%s - after LOL normalization", cond)
	plot(x.norm, main=str, xlim=Xlim, ylim=Ylim, xlab=x.lab, ylab="FC (log2)", pch=20, cex=0.7, col="grey")
	lines(c(0.8*Xlim[1], 1.2*Xlim[2]), c(0,0), col=4, lwd=1)
	lines(LWS.curve$x[order(LWS.curve$x)], LWS.curve$fitted[order(LWS.curve$x)], col=2, lwd=2)

        dev.off()

	x.norm

}

#-----------------------------------------------------#