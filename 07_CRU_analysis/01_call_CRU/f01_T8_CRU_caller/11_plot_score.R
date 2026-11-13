#!/usr/bin/env Rscript
##*## Score distribution for the scored cohort; descriptive figure only
suppressPackageStartupMessages({library(ggplot2); library(grid)})
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k,d=NULL) {i<-match(k,args);if(is.na(i))d else args[i+1L]}
scores <- read.delim(arg('--scores'),check.names=FALSE)
ref_path <- arg('--reference-scores')
out <- arg('--outdir');font <- arg('--font-family','Liberation Sans')
label <- arg('--label','CRU')
score <- scores$regulatory_score
stopifnot(nrow(scores)>0, all(is.finite(score)), all(score>=0&score<=100))
matched<-system2('fc-match',c('-f',shQuote('%{family}'),shQuote(font)),stdout=TRUE)
if(!any(grepl(tolower(font),tolower(matched),fixed=TRUE)))stop('Font unavailable: ',font)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
style<-theme_classic(base_size=11.5,base_family=font)+theme(
  panel.grid=element_blank(),text=element_text(colour='#222222'),
  axis.line=element_line(linewidth=.8),axis.ticks=element_line(linewidth=.8),
  plot.title=element_text(size=11.5,face='plain',hjust=.5),
  plot.tag=element_text(size=14,face='bold'),plot.tag.position=c(.01,.99),
  plot.margin=margin(8,9,6,6),legend.position='none')
a<-ggplot(data.frame(score),aes(score))+geom_histogram(binwidth=5,boundary=0,
  fill='#96CEDD',colour='white',linewidth=.35)+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=expansion(mult=c(.01,.01)))+
  scale_y_continuous(expand=expansion(mult=c(0,.06)))+
  labs(x='Regulatory score',y='CRU count',
       title=sprintf('%s score distribution (n = %s)',label,format(nrow(scores),big.mark=',')),
       tag='a')+style
panels<-list(a)
if(!is.null(ref_path) && file.exists(ref_path)){
  ref<-read.delim(ref_path,check.names=FALSE)
  d<-rbind(data.frame(dataset='Reference',score=ref$regulatory_score),
           data.frame(dataset=label,score=score))
  d$dataset<-factor(d$dataset,levels=c('Reference',label))
  labels<-data.frame(dataset=factor(c('Reference',label),levels=levels(d$dataset)),
                     label=c(sprintf('n = %s',format(nrow(ref),big.mark=',')),
                             sprintf('n = %s',format(nrow(scores),big.mark=','))))
  panels<-c(panels,list(ggplot(d,aes(dataset,score,fill=dataset))+geom_boxplot(width=.43,linewidth=.7,
    outlier.size=.5,outlier.alpha=.2)+scale_fill_manual(values=setNames(c('#96CEDD','#E89B6C'),levels(d$dataset)))+
    geom_text(data=labels,aes(x=dataset,y=104,label=label),inherit.aes=FALSE,family=font,size=3.3)+
    scale_y_continuous(limits=c(0,109),breaks=seq(0,100,25),expand=expansion(mult=c(0,0)))+
    labs(x=NULL,y='Regulatory score',title='Shared score scale',tag='b')+style))
}
pdf_file<-file.path(out,'01_score_distribution.pdf')
cairo_pdf(pdf_file,width=if(length(panels)>1) 6.8 else 4.2,height=3.5,family=font)
tryCatch({
  grobs<-lapply(panels,ggplotGrob)
  widths<-do.call(unit.pmax,lapply(grobs,function(g)g$widths))
  heights<-do.call(unit.pmax,lapply(grobs,function(g)g$heights))
  grid.newpage();pushViewport(viewport(layout=grid.layout(1,length(grobs))))
  for(i in seq_along(grobs)){grobs[[i]]$widths<-widths;grobs[[i]]$heights<-heights
    pushViewport(viewport(layout.pos.row=1,layout.pos.col=i));grid.draw(grobs[[i]]);popViewport()}
  popViewport()
},finally=dev.off())
cat('PDF:',pdf_file,'\n')
