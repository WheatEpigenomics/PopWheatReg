#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 03_plot_results.R all_chromosome_data.tsv pooled_data.tsv output_dir")
}

all_path <- args[[1]]
pooled_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

all_data <- read.delim(all_path, check.names=FALSE, stringsAsFactors=FALSE,
                       na.strings=c("NA","NaN","nan",""))
pooled_data <- read.delim(pooled_path, check.names=FALSE, stringsAsFactors=FALSE,
                          na.strings=c("NA","NaN","nan",""))
chromosomes <- unique(all_data$chrom)
pooled_label <- paste0("ALL_",length(chromosomes))
generations <- c(20,50,100,200)
blue <- "#2166AC"
red <- "#B2182B"
distance_bounds <- c(0,10000,20000,50000,100000,250000,500000,750000,1000000)
plot_audit <- list()

draw_panel <- function(dat, generation, chromosome) {
  z <- dat[dat$generation == generation,]
  z <- z[order(z$distance_low_bp),]
  if (nrow(z) != 8L || !identical(as.numeric(z$distance_low_bp),head(distance_bounds,-1)) ||
      !identical(as.numeric(z$distance_high_bp),tail(distance_bounds,-1))) {
    stop(paste("Expected all eight distance strata:",chromosome,generation))
  }
  valid <- z$n_supported_real_aglinks > 0
  numeric_fields <- c("real_mean_percent","real_sd_percent","real_q025_percent","real_q975_percent",
                      "random_mean_percent","random_sd_percent","random_q025_percent","random_q975_percent")
  if (any(!is.finite(as.matrix(z[valid,numeric_fields,drop=FALSE]))) ||
      any(z$real_sd_percent[valid] < 0) || any(z$random_sd_percent[valid] < 0)) {
    stop(paste("Unexpected non-finite value or negative SD:",chromosome,generation))
  }
  if (any(!is.na(as.matrix(z[!valid,numeric_fields,drop=FALSE])))) {
    stop(paste("Unsupported strata must have explicit NA values:",chromosome,generation))
  }
  # Distance strata are equally spaced categorical positions so both short
  # distance strata (0--10 and 10--20 kb) remain legible.
  x <- seq_len(8L)
  x_labels <- paste0(z$distance_low_bp/1000,"-",z$distance_high_bp/1000)
  lower <- c(z$real_mean_percent[valid]-z$real_sd_percent[valid],
             z$random_mean_percent[valid]-z$random_sd_percent[valid])
  upper <- c(z$real_mean_percent[valid]+z$real_sd_percent[valid],
             z$random_mean_percent[valid]+z$random_sd_percent[valid])
  min_value <- if (length(lower)) min(lower) else NA_real_
  max_value <- if (length(upper)) max(upper) else NA_real_
  use_50 <- generation %in% c(100,200) && length(lower)>0 && min_value>=50
  ymin <- if (use_50) 50 else if (length(lower)) min(0,floor(min_value/10)*10) else 0
  ymax <- if (length(upper)) max(105,ceiling(max_value/5)*5) else 105
  y_limits <- c(ymin,ymax)
  y_ticks <- if (use_50) seq(50,floor(ymax/10)*10,10) else pretty(y_limits,n=6)
  y_ticks <- y_ticks[y_ticks>=ymin & y_ticks<=ymax]
  plot(NA, xlim=c(.5,nrow(z)+.5), ylim=y_limits, axes=FALSE, xlab="", ylab="",xaxs="i")
  draw_sd <- function(means,sds,color) {
    lo <- means[valid]-sds[valid]; hi <- means[valid]+sds[valid]
    segments(x[valid],lo,x[valid],hi,col=color,lwd=.8)
    segments(x[valid]-.045,lo,x[valid]+.045,lo,col=color,lwd=.8)
    segments(x[valid]-.045,hi,x[valid]+.045,hi,col=color,lwd=.8)
  }
  draw_sd(z$real_mean_percent,z$real_sd_percent,blue)
  draw_sd(z$random_mean_percent,z$random_sd_percent,red)
  # NA rows retain their categorical position and break both connecting lines.
  lines(x,z$real_mean_percent,col=blue,lwd=1.5)
  points(x,z$real_mean_percent,col=blue,pch=16,cex=.75)
  lines(x,z$random_mean_percent,col=red,lwd=1.5)
  points(x,z$random_mean_percent,col=red,pch=16,cex=.75)
  axis(1,at=x,labels=FALSE)
  text(x,par("usr")[3]-0.035*diff(par("usr")[3:4]),labels=x_labels,
       srt=55,adj=1,xpd=TRUE,cex=1)
  axis(2,at=y_ticks,las=1,cex.axis=1)
  box(lwd=.7)
  title(main=paste("Generation",generation),cex.main=1,font.main=2,line=.7)
  if (!any(valid)) {
    label <- if (sum(z$n_all_real_aglinks)>0) "No pairs met the matching criteria" else "No eligible agLinks in these strata"
    text(4.5,mean(y_limits),label,cex=1)
  }
  plot_audit[[length(plot_audit)+1L]] <<- data.frame(
    chrom=chromosome,generation=generation,n_supported_strata=sum(valid),
    n_unsupported_strata=sum(!valid),unsupported_strata=paste(x_labels[!valid],collapse=","),
    y_min=ymin,y_max=ymax,min_mean_minus_sd=min_value,max_mean_plus_sd=max_value,
    retained_50_percent_lower_bound=use_50,error_bars_clipped=FALSE,
    stringsAsFactors=FALSE)
}

draw_figure <- function(dat, main_title, output_path, chromosome) {
  cairo_pdf(output_path,width=7.4,height=6.4,family="Arial",pointsize=6)
  on.exit(dev.off(),add=TRUE)
  layout(matrix(c(1,1,2,3,4,5),nrow=3,byrow=TRUE),heights=c(.24,1,1))
  par(oma=c(3.3,3.8,.3,.3),family="Arial",cex=1)
  par(mar=c(0,0,0,0))
  plot.new()
  text(.5,.80,main_title,cex=1,font=1)
  legend("center",legend=c("Real agLinks","Distance- and background-matched random pairs"),
         col=c(blue,red),lty=1,lwd=1.5,pch=16,pt.cex=.75,horiz=TRUE,
         bty="n",cex=1)
  for (generation in generations) {
    par(mar=c(4.2,3.4,2.2,.8))
    draw_panel(dat,generation,chromosome)
  }
  mtext("ACR-gene distance (kb)",side=1,outer=TRUE,line=1.5,cex=1)
  mtext("Within-LD-block links (%)",side=2,outer=TRUE,line=2.2,cex=1)
}

for (chromosome in chromosomes) {
  dat <- all_data[all_data$chrom == chromosome,]
  if (nrow(dat) == 0) stop(paste("No plotting data for chromosome",chromosome))
  draw_figure(dat,
              paste0("Chromosome ",chromosome,": real vs random agLinks across distance"),
              file.path(output_dir,paste0("Chromosome_",chromosome,"_distance_generations.pdf")),chromosome)
}

draw_figure(pooled_data,
            paste0("Pooled chromosomes (n=",length(chromosomes),"): real vs random agLinks across distance"),
            file.path(output_dir,paste0("All_",length(chromosomes),"_chromosomes_pooled_distance_generations.pdf")),pooled_label)

write.table(do.call(rbind,plot_audit),file.path(output_dir,"plot_limits_audit.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE,na="NA")

cat(sprintf("Created %d chromosome figures and one pooled figure in %s\n",length(chromosomes),output_dir))
