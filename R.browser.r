##########################################################################################
## Visualization
##########################################################################################

## Modules
##########################################################################################
data.loading <- function(datLink=NULL, dirName=NA, chrom=NA, start=NA, end=NA, labels=NA, nbins=end-start+1, name=NA, pattern="bdg", header=TRUE, sep="\t") {

  label.checking <- function(fileName, labels) {
    is.lable <- FALSE
    for(l in labels) {
      if(any(grepl(l, fileName, ignore.case=TRUE))) return(TRUE)
    }
    is.lable
  }

  if(is.null(datLink)){
    files <- dir(path=dirName, pattern=pattern, all.files=T, full.names=T)
    datList <- list(); datNames <- c()
    for(f in files) if(file.exists(f)) {
      dat <- read.table(f, header=header, sep=sep) # format: chrom[sep]start[sep]end[sep]value[sep]...
      datList <- c(datList, list(dat))
      datNames <- c(datNames, f)
    }
    names(datList) <- datNames
    save(list = c("datList", "files"), file = sprintf("%s/%s_tracking.RData", dirName, name))

    if(file.exists(sprintf("%s/%s_tracking.RData", dirName, name))) TRUE else FALSE

  } else {
    load(datLink)

    bins <- seq(start, end, length=nbins)
    visual.dat <- c(); lab.order <- c()
    for(f in files) {
      if(!label.checking(f, labels)) next

      dat <- datList[[f]]
      dat <- dat[dat[,1] == chrom & dat[,2] >= start & dat[,3] <= end,]
      dat <- dat[order(dat[,2]),]

      density <- rep(0, length(bins)-1); i <- c(1)
      for(j in 2:length(bins)) {
        while(i<=nrow(dat) && dat[i,3]<bins[j-1]) {
          i <- i + 1
          if(i>nrow(dat)) break
        }

        if(i>nrow(dat) || dat[i,2]>bins[j]) {
          density[j-1] <- c(-1)
          next
        }

        if(!(dat[i,3]<bins[j-1] || dat[i,2]>bins[j])) i.pre <- i

        tmp <- c()
        while(i<=nrow(dat) && !(dat[i,3]<bins[j-1] || dat[i,2]>bins[j])) {
          tmp <- c(tmp, dat[i,4])
          i <- i + 1
        }

        if(!is.null(tmp)) {
          density[j-1] <- mean(tmp)
          i <- i.pre
        }
      }
      visual.dat <- cbind(visual.dat, density)
      lab.order <- c(lab.order, f)
    }
    idx.lab <- 1:ncol(visual.dat); ord <- c(); labs <- c()
    for(i in 1:length(labels)) if(any(grepl(labels[i], lab.order))) {
       labs <- c(labs, labels[i])
       ord <- c(ord, idx.lab[grepl(labels[i], lab.order)])
    }

    list(dat=visual.dat, bin=bins, labOrder=ord, label=labs, chrom=chrom, start=start, end=end)
  }
}

browser <- function(visual.dat, bins, labels, colors, chrom, start, end, show.ylim=TRUE, layout.config=NA, is.smooth=FALSE, is.commmon=FALSE, ylim=c(0, max(visual.dat)*1.3), nbins=1000, cex=1, main=NULL,
  xlab=NULL, at=NULL, layout.mode="separate", scale.show=5000, species="mm10", layout.width=c(1,10), scale.signal=1) {
      
  chrom.id <- sprintf("%s:", chrom)
  if(scale.signal != 1) {
    visual.dat[visual.dat>0] <- visual.dat[visual.dat>0] / scale.signal
  }
  
  if(!is.matrix(ylim) && length(ylim)<=2) ylim <- matrix(c(rep(0, ncol(visual.dat)), apply(visual.dat, 2, max, na.rm=T)), ncol=2) # matrix(rep(ylim, ncol(visual.dat)), ncol=2, byrow=T)
  if(length(show.ylim)!=ncol(visual.dat)) show.ylim <- rep(TRUE, ncol(visual.dat))
  if(length(is.smooth)!=ncol(visual.dat)) is.smooth <- rep(FALSE, ncol(visual.dat))
  
  x.axis <- (bins[2:length(bins)] + bins[1:(length(bins)-1)]) / 2
  if(is.commmon) {
    x.axis <- x.axis[rowSums(visual.dat<0)<=0]
    visual.dat <- visual.dat[rowSums(visual.dat<0)<=0,]
  }

  if(any(is.na(layout.config))) {layout(cbind(rep(0,ncol(visual.dat)+1), matrix(1:(ncol(visual.dat)+1),ncol=1)), w=layout.width, h=c(5, rep(3,ncol(visual.dat))))} else {layout(cbind(rep(0,ncol(visual.dat)+1), matrix(1:(ncol(visual.dat)+1),ncol=1)), w=layout.width, h=c(5, layout.config))}
  
  par(mar=c(1,1,6,10)*0.5)
  plot(NA, bg='gray', ann=F, ylim=c(0, 1), xlim=c(start, end), cex=cex, cex.lab=cex, cex.main=cex, xaxs='i', axes=F)
  lines(seq(start+0.3*(end-start), start+0.3*(end-start)+scale.show, len=5), rep(0.5,5), lwd=cex/2)
  lines(rep(start+0.3*(end-start),3), c(0, 0.1, 0.2)+0.4, lwd=cex/2)
  lines(rep(start+0.3*(end-start)+scale.show,3), c(0, 0.1, 0.2)+0.4, lwd=cex/2)
  if(scale.show > 1000) text(start+0.3*(end-start), 0.5, labels=sprintf("%.fkb", scale.show/1000), cex=0.5*cex, pos=2) else text(start+0.3*(end-start), 0.5, labels=sprintf("%.fbp", scale.show), cex=0.5*cex, pos=2)
  text(start+0.3*(end-start)+scale.show, 0.5, labels=species, cex=0.5*cex, pos=4)
  axis(3, at=seq(start,end,len=5), label=print.money(seq(start,end,len=5)), cex.axis=0.5*cex, las=1, lwd=0, lwd.ticks=cex/2)
  axis(2, at=par()$usr[4], labels=chrom.id, las=2, cex.axis=0.5*cex, lwd=2, font.axis=2, col.axis="black", tick=F, hadj=1+min(c(nchar(chrom.id), nchar(labels)))/nchar(chrom.id))
  
  for(i in 1:ncol(visual.dat)) {
    density <- visual.dat[,i]
    if(is.smooth[i]) fit <- smooth.spline(x.axis, density, df=ifelse(length(x.axis) > nbins, nbins, length(x.axis))) else fit <- list(x=x.axis, y=density)
    fit$y[fit$y<0] <- 0
    if(all(ylim[i,] == 0)) ylim[i,] <- c(0, max(visual.dat)*1.3)
    
    if(layout.mode == "collapse") {
      if(i == 1) {
        par(mar=c(0,1,1,10)*0.5)
        plot(NA, bg='gray', ann=F, yaxt='n', xaxt='n', ylim=ylim[i,], xlim=c(start, end), frame.plot=T, cex=cex, cex.lab=cex, cex.main=cex, xaxs='i')
      } else if(i == ncol(visual.dat)){
        par(mar=c(1,1,0,10)*0.5)
        plot(NA, bg='gray', ann=F, yaxt='n', xaxt='n', ylim=ylim[i,], xlim=c(start, end), frame.plot=T, cex=cex, cex.lab=cex, cex.main=cex, xaxs='i')
      } else {
        par(mar=c(0,1,0,10)*0.5)
        plot(NA, bg='gray', ann=F, yaxt='n', xaxt='n', ylim=ylim[i,], xlim=c(start, end), frame.plot=T, cex=cex, cex.lab=cex, cex.main=cex, xaxs='i')
      }
      axis(2, at=ylim[i,2]/2, labels=labels[i], las=2, cex.lab=cex, cex.axis=cex, lwd=2, font.axis=2, col.axis=colors[i])
      if(show.ylim[i]) text(start+0.05*(end-start), 0.9*ylim[i,2], ifelse(ylim[i,2]>1000 || ylim[i,2] < 1e-2, format(ylim[i,2], sci=T, digits=2), sprintf("[0,%.2f]", ylim[i,2])), cex=0.8*cex, font=2)
      if(i == ncol(visual.dat) && !is.null(xlab)) axis(1, at=at, labels=xlab, cex.axis=cex, cex.lab=cex, cex=cex)
      if(i == 1 && !is.null(main)) axis(3, at=start+(end-start)/2, tick=F, labels=sprintf("%s (%s:%s-%s)", main, chrom, start, end), cex.axis=cex, cex.lab=cex, cex=cex)
    } else if(layout.mode == "separate") {
        par(mar=c(1,1,1,10)*0.5)
        plot(NA, bg='gray', ann=F, ylim=ylim[i,], xlim=c(start, end), cex=cex, cex.lab=cex, cex.main=cex, xaxs='i', axes=F)
        if(show.ylim[i]) axis(2, at=ylim[i,], labels=c("0", ifelse(ylim[i,2]>1000 || ylim[i,2] < 1e-2, format(ylim[i,2], sci=T, digits=2), sprintf("%.2f", ylim[i,2]))), cex.axis=0.5*cex, las=2, col.axis=colors[i], col=colors[i])
        axis(2, at=ylim[i,2]/2, labels=labels[i], las=2, cex.axis=0.5*cex, lwd=2, font.axis=2, col.axis=colors[i], col=colors[i], tick=F, hadj=1+min(c(nchar(chrom.id), nchar(labels)))/nchar(labels[i]))
    }
    
    if(is.smooth[i]) polygon(c(rev(fit$x), fit$x), c(rev(fit$y), 0*fit$y), col = colors[i], border = NA) else {
      for(k in 1:length(fit$x)) {
        if(!show.ylim[i] && fit$y[k] <= 0) col <- "transparent" else col <- colors[i]
        if(fit$y[k]>=0) {
          lines(c(fit$x[k],fit$x[k]), c(0, fit$y[k]), col = col, lwd=1)
        } else {
          lines(c(fit$x[k],fit$x[k]), c(fit$y[k], 0), col = col, lwd=1)
        }
      }
    }
  }
  
}

.screener <- function(dirName, chrom, start, end, tmpDir=sprintf("%s/tmp_%.4f", dirName, runif(1)), pattern="bdg") {
  if(!dir.exists(tmpDir)) dir.create(tmpDir)

  for(f in dir(path=dirName, pattern=pattern, full.names=T)) {
    outFile <- sprintf("%s/%s", tmpDir, gsub(".*/", "", f))
    if(grepl("gz$", f)) {
      system(sprintf("gunzip -c %s | awk '{if($1==\"%s\" && $2>=%d && $3<=%d) printf(\"%%s\\t%%s\\t%%s\\t%%s\\n\", $1, $2, $3, $4);}' | gzip -c > %s", f, chrom, start, end, outFile))
    } else {
      system(sprintf("cat %s | awk '{if($1==\"%s\" && $2>=%d && $3<=%d) printf(\"%%s\\t%%s\\t%%s\\t%%s\\n\", $1, $2, $3, $4);}' > %s", f, chrom, start, end, outFile))
    }
  }

  return(tmpDir)
}

.plot.visual <- function(dirName="data", tmpDir=NA, outdir=NA, show.ylim=FALSE, layout.config=NA, colors=NA, chrom=NA, start=NA, end=NA, is.commmon=FALSE, is.bigDat=FALSE, layout.mode="separate", scale.show=5000, cex=3.5, species="mm10",
  pattern="bdg.gz", labels=NA, header=FALSE, name=NA, is.smooth=FALSE, nbins=1000, main=NULL, xlab=NULL, at=NULL, fwidth=2000, fheight=500, layout.width=c(1,10), scale.signal=1) {
  ## import data
  if(is.bigDat && !file.exists(sprintf("%s/%s_tracking.RData", dirName, name))) {
    tmpDir <- .screener(dirName, chrom, start, end, pattern=pattern)
    if(is.na(outdir)) outdir <- dirName
    dirName <- tmpDir
  }

  datLink <- sprintf("%s/%s_tracking.RData", dirName, name)
  if(!file.exists(datLink)) {
    data.loading(dirName=dirName, header=FALSE, pattern=pattern, name=name)
  }
  mat <- data.loading(datLink=datLink, chrom=chrom, start=start, end=end, labels=labels, header=FALSE, name=name)

  ## plot
  if(is.na(outdir)) bmp(sprintf("%s/%s.tracking.tiff", dirName, name), width = fwidth, height = fheight) else bmp(sprintf("%s/%s.tracking.tiff", outdir, name), width = fwidth, height = fheight)
  if(all(is.na(layout.config))) layout.config <- rep(3,ncol(mat$dat))
  browser(visual.dat=mat$dat[,mat$labOrder], bins=mat$bin, show.ylim=show.ylim, layout.config=layout.config, labels=labels, is.smooth=is.smooth, is.commmon=is.commmon, colors=colors, chrom=chrom, start=start,
    layout.mode=layout.mode, end=end, scale.show=scale.show, cex=cex, species=species, layout.width=layout.width, scale.signal=scale.signal)
  dev.off()
  
  if(file.exists(sprintf("%s/%s_tracking.RData", tmpDir, name))) system(sprintf("mv %s/%s_tracking.RData %s", tmpDir, name, outdir))
  if(dir.exists(sprintf("%s", tmpDir))) unlink(tmpDir, recursive=TRUE, force=TRUE)
}

makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

print.money <- function(x, ...) {
  paste0(formatC(as.double(x), format="f", digits=0, big.mark=","))
}

makeColors <- function(files, labels, sep=".*/|rep.*", alpha=1) {
  colbars <- gsub(sep, "", files)
  warning(sprintf("The default colors will be grouped by labels: %s", paste0(colbars, collapse=",")))

  labels <- labels[order(colbars)]
  colbars <- colbars[order(colbars)]
  cols <- rainbow(length(unique(colbars)))
  names(cols) <- unique(colbars)
  colors <- colbars
  for(lab in unique(colbars)) {
    colors[colbars %in% lab] <- makeTransparent(cols[lab], alpha=alpha)
  }
    
  colors
}

makeScale_Show <- function(n) {
  scale <- as.integer(n/4)
  if(scale < 1000 && scale >= 100) {
    round(scale / 100) * 100
  } else if(scale >= 1000 && scale < 10000) {
    round(scale / 1000) * 1000
  } else if(scale >= 10000) {
    round(scale / 5000) * 5000
  } else {
    50
  }
}

makeLabelsColors <- function(opt, name) {
  if(is.na(opt$labels)) {
    labels <- gsub(sprintf(".*/|%s.*", opt$pattern), "", dir(opt$directory, pattern=opt$pattern, full.names=T))
  } else {
    labels <- unlist(strsplit(opt$labels, ","))
    label.check <- rep(FALSE, length(labels))
    names(label.check) <- labels
    Files <- dir(opt$directory, pattern=opt$pattern, full.names=T)
    for(lab in labels) {
      if(length(Files) == length(labels)) {
        if(any(grepl(lab, Files))) label.check[lab] <- TRUE
      } else if(file.exists(sprintf("%s/%s_tracking.RData", opt$directory, name))) {
        load(sprintf("%s/%s_tracking.RData", opt$directory, name))
        if(any(grepl(lab, files))) label.check[lab] <- TRUE
      }
    }
    if(!all(label.check)) {
      warning(sprintf("Invalid labels: %s", opt$labels))
      labels <- gsub(sprintf(".*/|%s.*", opt$pattern), "", dir(opt$directory, pattern=opt$pattern, full.names=T))
      warning(sprintf("The default labels will be used: %s", paste0(labels, collapse=",")))
    }
  }

  if(is.na(opt$colors)) {
    colors <- makeColors(dir(opt$directory, pattern=opt$pattern, full.names=T), labels, sep=opt$colorPattern, alpha=opt$color_transparent)
  } else {
    colors <- unlist(strsplit(opt$colors, ","))
    if(length(colors) != length(labels)) {
      colors <- makeColors(dir(opt$directory, pattern=opt$pattern, full.names=T), labels, sep=opt$colorPattern, alpha=opt$color_transparent)
    }
  }

  if(length(unique(labels)) != length(labels)){
    stop("Must have the unique labels for all tracks")
  }
  
  if(length(colors) != length(labels)){
    warning("Number of colors different from the number of labels")
    if(length(colors) > length(labels)) colors <- colors[1:length(labels)]
    if(length(colors) < length(labels)) colors <- rep(colors, ceiling(length(labels)/length(colors)))[1:length(labels)]
  }
  
  list(labels=labels, colors=colors)
}
##########################################################################################

###############################
local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.cnr.berkeley.edu/"
  options(repos = r)
})

tryCatch(
  library("optparse"),
  error=function(e) {
    print(e)
    install.packages("optparse", lib=.libPaths()[1L], quiet=F)
  },
  finally=library("optparse")
)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, help="To specify a directory which stores BEDGRAPH files. Default: NULL. Required.", metavar="character"),
  make_option(c("-o", "--out_directory"), type="character", default=NA, help="To specify an output directory. Default: the same with directory.", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="single", help="To specify a way for single task or a batch of tasks (single or batch). Default: single.", metavar="character"),
  make_option(c("-f", "--filename"), type="character", default=NULL, help="To specify a filenames to store a series of genomic coordinates, such as: 'chrom start end...' by TAB. Default: NULL.", metavar="character"),
  make_option(c("-c", "--chrom"), type="character", default=NULL, help="To specify a chromosome name. Default: NULL.", metavar="character"),
  make_option(c("-s", "--start"), type="integer", default=NULL, help="To specify a start of genomic coordinate. Default: NULL.", metavar="character"),
  make_option(c("-e", "--end"), type="integer", default=NULL, help="To specify an end of genomic coordinate. Default: NULL.", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="bdg", help="To specify a string pattern embedding in file names of BEDGRAPH files. Default: bdg.", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NA, help="To specify a string pattern in the output figures. Default: Region.chrom_start_end.", metavar="character"),
  make_option(c("--show_ylim"), type="character", default=NA, help="To specify whether to show ylim. Default: TRUE for each.", metavar="character"),
  make_option(c("--is_commmon_ylim"), type="logical", default=FALSE, help="To specify whether to draw the same scale of ylim. Default: FALSE.", metavar="character"),
  make_option(c("--is_bigData"), type="logical", default=TRUE, help="To specify whether to use AWK mode. Default: TRUE.", metavar="character"),
  make_option(c("--layout_mode"), type="character", default="separate", help="To specify which layout mode will be used (speparate or collapse). Default: speparate.", metavar="character"),
  make_option(c("--scale_show"), type="integer", default=5000, help="To specify the length to show in visualization figure. Default: 5000 (that is 5kb).", metavar="character"),
  make_option(c("--cex"), type="integer", default=3.5, help="To specify the size to show labels. Default: 3.5.", metavar="character"),
  make_option(c("--figure_width"), type="integer", default=2500, help="To specify the width of figure. Default: 2500.", metavar="character"),
  make_option(c("--figure_height"), type="integer", default=600, help="To specify the height of figure. Default: 600.", metavar="character"),
  make_option(c("--is_smooth"), type="logical", default=FALSE, help="To specify whether to smooth read signals. Default: FALSE.", metavar="character"),
  make_option(c("--labels"), type="character", default=NA, help="To specify the labels to show; format: lab1,lab2,...labN. Default: the labels are splitted by 'PATTEN'.", metavar="character"),
  make_option(c("--colors"), type="character", default=NA, help="To specify the colors corresponding to labels; format: color1,color2,...colorN. Default: colors will be grouped by labels (split by 'rep').", metavar="character"),
  make_option(c("--species"), type="character", default="mm10", help="To specify the species name to show. Default: mm10.", metavar="character"),
  make_option(c("--colorPattern"), type="character", default=".*/|rep.*", help="To specify a pattern to split file names to group them for showing. Default:'.*/|rep.*'.", metavar="character"),
  make_option(c("--color_transparent"), type="double", default=1, help="To specify the degree of transparent color between 0 (transparent) and 1. Default: 1.", metavar="character"),
  make_option(c("--layout_width"), type="character", default="1:10", help="To specify the ratio of the label panel over the figure panel. Default: '1:10'.", metavar="character"),
  make_option(c("--scale_signal"), type="integer", default=1, help="To specify an integer factor to scale the read signal data. Default: '1:10'.", metavar="character"),
  make_option(c("--layout_config"), type="character", default=NA, help="To specify the height for each track. Default: 3 for each.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
### configuration
if(is.null(opt$directory)) {
  stop("Must specifiy a directory which stores BEDGRAPH files.")
}

dirName <- opt$directory
trim_pattern <- opt$pattern
pattern <- opt$pattern

if(is.na(opt$out_directory)) {
  outdir <- dirName
} else {
  outdir <- opt$out_directory
  if(!dir.exists(outdir)) dir.create(outdir)
}

if(is.na(opt$show_ylim)) show.ylim <- TRUE else {
  show.ylim <- as.logical(unlist(strsplit(opt$show_ylim, ",")))
}
is.bigDat <- as.logical(opt$is_bigData)
if(is.bigDat) warning("To extract data from BEDGRAPH files based on AWK environment")
is.commmon <- as.logical(opt$is_commmon_ylim)
cex <- opt$cex
scale.show <- opt$scale_show
fwidth <- opt$figure_width
fheight <- opt$figure_height
is.smooth <- as.logical(opt$is_smooth)
layout.mode <- opt$layout_mode
species <- opt$species
layout.width <- as.double(unlist(strsplit(opt$layout_width, ":")))
scale.signal <- as.double(opt$scale_signal)

if(opt$mode == "single") {
  if(is.null(opt$chrom) || is.null(opt$chrom) || is.null(opt$end)) stop ("Must specify a valid genomic coordinate")
  chrom <- opt$chrom
  start <- as.double(opt$start)
  end <- as.double(opt$end)
  
  if(!is.na(opt$name)) {
    name <- opt$name
  } else {
    name <- sprintf("Region.%s_%s_%s", chrom, start, end)
  }
  
  labCol.list <- makeLabelsColors(opt, name)
  labels <- labCol.list$labels
  colors <- labCol.list$colors
  if(is.na(opt$layout_config)) layout.config <- rep(3,length(labels)) else {
    layout.config <- as.integer(unlist(strsplit(opt$layout_config, ",")))
  }
  
  # to draw figure
  .plot.visual(dirName=dirName, outdir=outdir, chrom=chrom, start=start, end=end, show.ylim=show.ylim, layout.config=layout.config, labels=labels, is.bigDat=is.bigDat, pattern=pattern, cex=cex, scale.show=scale.show,
    name=name, colors=colors, is.commmon=is.commmon, fwidth=fwidth, fheight=fheight, is.smooth=is.smooth, layout.mode=layout.mode, species=species, layout.width=layout.width, scale.signal=scale.signal)
} else if(opt$mode == "batch") {
  if(!file.exists(opt$filename)) stop(sprintf("Invalid file name: %s", opt$filename))
  warning("By default, there is a file header: chrom start end ...")
  regions <- read.table(opt$filename, header=T)
  if(ncol(regions)>=4) warning("By default, the 4th column will be uased as a part of the output file name.")
  
  for(i in 1:nrow(regions)) {
    chrom <- paste0(regions[i,1])
    start <- as.double(regions[i,2])
    end <- as.double(regions[i,3])
    if(ncol(regions)>=4) {
      name <- paste0(regions[i,4])
    } else {
      name <- sprintf("Region.%s_%s_%s", chrom, start, end)
    }
    
    var.logical <- file.exists(sprintf("%s/%s_tracking.RData", dirName, name)) ||
      length(dir(dirName, pattern=pattern, full.names=T)) > 0

    # to draw figure
    if(var.logical) {
      labCol.list <- makeLabelsColors(opt, name)
      labels <- labCol.list$labels
      colors <- labCol.list$colors
      if(is.na(opt$layout_config)) layout.config <- rep(3,length(labels)) else {
        layout.config <- as.integer(unlist(strsplit(opt$layout_config, ",")))
      }
    
      scale.show <- makeScale_Show(end - start)
      warning(sprintf("The reset scale showing in figure when using a BATCH mode: %s", scale.show))
      
      .plot.visual(dirName=dirName, outdir=outdir, chrom=chrom, start=start, end=end, show.ylim=show.ylim, layout.config=layout.config, labels=labels, is.bigDat=is.bigDat, pattern=pattern, cex=cex, scale.show=scale.show,
        name=name, colors=colors, is.commmon=is.commmon, fwidth=fwidth, fheight=fheight, is.smooth=is.smooth, layout.mode=layout.mode, species=species, layout.width=layout.width, scale.signal=scale.signal)
    } else {
      next
    }
  }
}

for(d in dir(path=dirName, pattern="tmp_", full.name=F)) {
  unlink(d, recursive=TRUE, force=TRUE)
}