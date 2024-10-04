args <- commandArgs(trailingOnly=T)
nucmer_show <- args[1] # nucmer show-coordi output
qryname <- args[2]
refname <- args[3]
qseqrm <- args[4]
rseqrm <- args[5]
linewidth <- args[6]
expandoff <- args[7] # no or yes
outpdf <- args[8] # PDF output file

linewidth <- as.numeric(linewidth)

if (expandoff == "yes") {
	lend.turnoff <- TRUE
} else {
	lend.turnoff <- FALSE
}

dotplot <- function(datafile, lend.turnoff=F,
                    refname="ref", qryname="asm",
					xlabel.rm = "chr", ylabel.rm = "tig0+",
                    line.width.factor=3.5,
					outpdf) {                  
  stopifnot(file.exists(datafile))	
  pdf.width = 6
  pdf.height = 6
  co <- read.delim(datafile, stringsAsFactor=F)
  if (nrow(co)<1) {
  	stop("no Nucmer alignments passing the criteria")
  }
  co <- read.delim(datafile, stringsAsFactor=F)
  print(head(co))
  colnames(co) <- c("firsts", "firste", "seconds", "seconde",
                    "first.align.l", "second.align.l", "ident",
                    "firstl", "secondl", "first.cov",
                    "second.cov", "first", "second")
  co$firsts <- as.numeric(as.character(co$firsts))
  co <- co[order(co$first, co$firsts, co$second, co$seconds), ]

  ### determine the total size of all the contigs in each set
  co$first <- as.character(co$first)
  co$second <- as.character(co$second)
  first.contigs.size <- tapply(co$firstl, co$first, max)
  second.contigs.size <- tapply(co$secondl, co$second, max)
  second.contigs.size <- second.contigs.size[unique(co$second)]
  
  print(second.contigs.size)
  
  ### maximum values for each contig:
  nfirst <- length(first.contigs.size)
  first.accum <- rep(0, times=nfirst)
  if (nfirst>1) {
    for (i in 2:nfirst) {
      first.accum[i] <- sum(first.contigs.size[1:(i-1)])   
    }
  } else {
    first.accum <- 0
  }
  names(first.accum) <- names(first.contigs.size)

  ### second set:
  nsecond <- length(second.contigs.size)
  second.accum <- rep(0, times=nsecond)
  if (nsecond>1) {
    for (i in 2:nsecond) {
      second.accum[i] <- sum(second.contigs.size[1:(i-1)])   
    }
  } else {
    second.accum <- 0
  }
  names(second.accum) <- names(second.contigs.size)
  
  max.size <- max(sum(first.contigs.size), sum(second.contigs.size))
  max.xsize <- max(sum(first.contigs.size))
  max.ysize <- max(sum(second.contigs.size))

  co$first.accum.s <- co$firsts + as.numeric(first.accum[co$first])
  co$first.accum.e <- co$firste + as.numeric(first.accum[co$first])
  co$second.accum.s <- co$seconds + as.numeric(second.accum[co$second])
  co$second.accum.e <- co$seconde + as.numeric(second.accum[co$second])

  ### plot
  pdf(outpdf, width=pdf.width, height=pdf.height)
  
  xrange <- c(-sum(first.contigs.size) / 50, sum(first.contigs.size))
  yrange <- c(0, sum(second.contigs.size))

  ### plot
  plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange,
       xlab=refname, ylab=qryname, cex.lab=1.2)
  col1 <- rgb(0, 0.4, 0, 0.5)
  col2 <- rgb(1, 0, 0, 0.5)
  ctg.num <- length(first.contigs.size)
  col.db <- data.frame(Col=rep(c(col1, col2), ctg.num)[1:ctg.num], Contig=names(first.contigs.size))
  
  abline(v=first.accum[1], lwd=0.6)
  abline(v=first.accum[-1], lwd=1, col="light grey")
  abline(v=sum(first.contigs.size), lwd=0.6)
  
  abline(h=second.accum[1], lwd=0.6)
  abline(h=second.accum, lwd=1, col="light grey")
  abline(h=sum(second.contigs.size), lwd=0.6)
  
  for (i in 1:nrow(co)) {
    plot.col <- as.character(col.db$Col[col.db$Contig==co[i, "first"]])
    distance <- co$second.accum.e[i] - co$second.accum.s[i]
    lend.val <- 2
    if (distance > max.size/50) {
      lend.val <- 1
    }
    if (lend.turnoff) {
      lend.val <- 1
    }
    
    lines(c(co$first.accum.s[i], co$first.accum.e[i]),
          c(co$second.accum.s[i], co$second.accum.e[i]),
          lwd=lend.val*line.width.factor, col=plot.col, lend=lend.val)
  }
  
  first.coord <- (c(first.accum[-1], sum(first.contigs.size)) + first.accum) / 2
  cat(first.coord)
  second.coord <- (c(second.accum[-1], sum(second.contigs.size)) + second.accum) / 2
  cat(second.coord)
  xlabels <- names(first.accum)
  xlabels <- gsub(xlabel.rm, "", xlabels)
  ylabels <- names(second.accum)
  cat("--", ylabels, "\n")
  ylabels <- gsub(ylabel.rm, "", ylabels)
  cat("==", ylabels, "\n")
  text(x=first.coord, y= - max.ysize / 50, labels=xlabels, cex = 0.8, xpd = T)
  text(x= - max.xsize / 40, y=second.coord, labels=ylabels, cex = 0.8, xpd = T)
  
  dev.off()
}

# plot
dotplot(datafile=nucmer_show, refname=refname, qryname=qryname,
        xlabel.rm=rseqrm, ylabel.rm=qseqrm,
		lend.turnoff=lend.turnoff, line.width.factor=linewidth,
		outpdf=outpdf)

