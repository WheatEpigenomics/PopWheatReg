#!/usr/bin/env Rscript
##*## One chromosome, saved calls only; fixed-window vector PDF and PNG atlas
args <- commandArgs(trailingOnly=TRUE)
options(scipen=999)
arg <- function(k,d=NULL) { i<-match(k,args); if(is.na(i)) d else args[i+1L] }
root<-arg('--chrom-dir'); inputs<-arg('--input-dir'); chrom<-arg('--chrom')
size<-as.numeric(arg('--size')); window<-as.numeric(arg('--window-bp','5000000'))
font<-arg('--font-family','Liberation Sans'); dpi<-as.integer(arg('--dpi','120'))
cap<-as.integer(arg('--max-aglink-arcs','1200')); seed<-as.integer(arg('--seed','20260805'))
co_cap<-as.integer(arg('--max-coacr-arcs','1200'))
stopifnot(!is.null(root),!is.null(inputs),!is.null(chrom),is.finite(size),size>0,window>0,dpi>0,cap>=0)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
family<-system2('fc-match',c('-f',shQuote('%{family}'),shQuote(font)),stdout=TRUE)
if(!any(grepl(tolower(font),tolower(family),fixed=TRUE))) stop('Requested font is unavailable: ',font)
tab<-function(p) {
  con<-if(grepl('\\.gz$',p))gzfile(p,'rt')else file(p,'rt')
  on.exit(close(con));read.delim(con,check.names=FALSE)
}
bed<-function(p) {
  if(file.info(p)$size==0)return(data.frame(chrom=character(),start=numeric(),end=numeric(),id=character()))
  d<-read.table(p,header=FALSE,comment.char='#');names(d)[1:3]<-c('chrom','start','end')
  if(ncol(d)>=4)names(d)[4]<-'id'
  stopifnot(all(d$chrom==chrom),all(d$start>=0),all(d$end<=size),all(d$end>d$start));d
}
density_path<-file.path(root,'f04_plot_data/01_density.tsv.gz')
edges_path<-file.path(root,'f04_plot_data/02_aglink.tsv.gz')
co_path<-file.path(inputs,'08_coacr.tsv')
cru_path<-file.path(root,'f02_calling',paste0(chrom,'.CRU.bed'))
block_path<-file.path(inputs,'05_blocks.bed')
paths<-c(density_path,edges_path,co_path,cru_path,block_path)
before<-tools::md5sum(paths)
density<-tab(density_path);es<-tab(edges_path);co<-tab(co_path)
crus<-bed(cru_path);blocks<-bed(block_path)
stopifnot(!any(crus$start[-1]<head(crus$end,-1)))
starts<-seq(0,size-1,by=window);ends<-pmin(size,starts+window)
stopifnot(starts[1]==0,tail(ends,1)==size,all(head(ends,-1)==tail(starts,-1)))
out<-file.path(root,'f05_figures');dir.create(out,showWarnings=FALSE)
png_dir<-file.path(out,'f01_png');dir.create(png_dir,showWarnings=FALSE)
atlas<-file.path(out,paste0(chrom,'.CRU.',format(window/1e6,trim=TRUE),'Mb.pdf'))
temporary<-paste0(atlas,'.tmp.pdf')
style<-theme_classic(base_size=12,base_family=font)+theme(
  panel.grid=element_blank(),text=element_text(colour='#222222'),axis.line=element_blank(),
  axis.ticks=element_line(linewidth=.8),panel.border=element_rect(fill=NA,colour='black',linewidth=.85),
  legend.position='top',legend.title=element_blank(),legend.text=element_text(size=10),
  plot.margin=margin(3,8,3,6),axis.title.y=element_text(size=10),plot.title=element_text(size=13))
hidex<-theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank())
interval_track<-function(d,a,b,label,colour,last=FALSE) {
  d<-d[d$start<b&d$end>a,,drop=FALSE]
  p<-ggplot(d)+geom_rect(aes(xmin=start/1e6,xmax=end/1e6,ymin=.25,ymax=.75),
    fill=colour,colour='#333333',linewidth=.7)+coord_cartesian(xlim=c(a,b)/1e6,ylim=c(0,1),expand=FALSE)+
    labs(x=if(last)paste0(chrom,' coordinate (Mb)')else NULL,y=label)+style+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  if(!last)p<-p+hidex
  p
}
index<-list()
arc_track<-function(edges,a,b,limit,label,colour,offset=0,heading=NULL) {
  use<-edges[edges$left<b&edges$right>=a,,drop=FALSE];original_n<-nrow(use)
  set.seed(as.integer((seed+a/5000+sum(utf8ToInt(chrom))*10000+offset)%%.Machine$integer.max))
  if(limit>0&&nrow(use)>limit)use<-use[sort(sample.int(nrow(use),limit)),,drop=FALSE]
  curves<-data.frame(x=numeric(),y=numeric(),edge=integer())
  if(nrow(use))curves<-do.call(rbind,lapply(seq_len(nrow(use)),function(j) {
    x<-seq(max(a,use$left[j]),min(b,use$right[j]),length.out=32)
    t<-(x-use$left[j])/max(1,use$right[j]-use$left[j])
    data.frame(x=x/1e6,y=sqrt(max(1,use$right[j]-use$left[j]))*4*t*(1-t),edge=j)
  }))
  top<-if(nrow(curves))max(1,max(curves$y))/.8 else 1
  p<-ggplot(curves,aes(x,y,group=edge))+
    geom_line(colour=colour,linewidth=.55,alpha=max(.04,min(.25,5/sqrt(max(1,nrow(use))))))+
    coord_cartesian(xlim=c(a,b)/1e6,ylim=c(0,top),expand=FALSE)+
    labs(x=NULL,y=paste0(label,'\nn=',original_n),title=heading)+style+hidex+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  list(plot=p,total=original_n,drawn=nrow(use))
}
##**# Draw on one persistent PDF device; PNG devices must not add PDF pages
cairo_pdf(temporary,width=14,height=8,family=font,onefile=TRUE)
pdf_device<-dev.cur()
tryCatch({for(i in seq_along(starts)) {
  a<-starts[i];b<-ends[i];xx<-c(a,b)/1e6
  ag<-arc_track(es,a,b,cap,'agLINK','#2E8DC8',heading=sprintf('%s:%g-%g Mb',chrom,a/1e6,b/1e6))
  ca<-arc_track(co,a,b,co_cap,'coACR','#9D83B9',offset=7919)
  panels<-list(ag$plot,ca$plot)
  dd<-density[density$position>=a-5000&density$position<=b+5000,,drop=FALSE]
  keys<-c('raw','density25','density105','density225','background')
  ll<-data.frame(x=rep(dd$position/1e6,length(keys)),y=unlist(dd[keys],use.names=FALSE),
    signal=factor(rep(keys,each=nrow(dd)),levels=keys))
  panels[[3]]<-ggplot(ll,aes(x,y,colour=signal))+geom_line(linewidth=.9,alpha=.85)+
    scale_colour_manual(values=c(raw='#C1C0BF',density25='#96CEDD',density105='#2E8DC8',density225='#5FB6BD',background='#DE6E68'),
      labels=c('Full raw span','25-kb smoothing','105-kb smoothing','225-kb smoothing','Local background'))+
    scale_y_continuous(expand=expansion(mult=c(.06,.12)))+coord_cartesian(xlim=xx)+
    labs(x=NULL,y='Span coverage')+style+hidex
  panels[[4]]<-interval_track(crus,a,b,'CRU','#E89B6C')
  panels[[5]]<-interval_track(blocks,a,b,'Haplotype blocks','#D3D3D3',TRUE)
  grobs<-lapply(panels,ggplotGrob);widths<-do.call(unit.pmax,lapply(grobs,function(g)g$widths))
  weights<-c(.85,.85,1.05,.45,.45)
  for(j in seq_along(grobs)) {
    grobs[[j]]$widths<-widths
    rr<-unique(grobs[[j]]$layout$t[grepl('^panel',grobs[[j]]$layout$name)])
    grobs[[j]]$heights[rr]<-unit(weights[j],'null')
  }
  g<-Reduce(function(x,y)rbind(x,y,size='first'),grobs)
  dev.set(pdf_device)
  if(i>1L)grid.newpage()
  grid.draw(g)
  preview<-file.path(png_dir,sprintf('%s_%d_%d.png',chrom,a,b))
  png(preview,width=14,height=8,units='in',res=dpi,type='cairo')
  tryCatch({grid.newpage();grid.draw(g)},finally=dev.off())
  index[[i]]<-data.frame(chrom=chrom,start=a,end=b,page=i,aglink_n=ag$total,drawn_n=ag$drawn,
    coacr_n=ca$total,drawn_coacr_n=ca$drawn,
    cru_n=sum(crus$start<b&crus$end>a),
    coverage_fraction=sum(pmax(0,pmin(crus$end,b)-pmax(crus$start,a)))/(b-a),pdf=atlas,png=preview)
  if(i%%10L==0L)cat(chrom,i,'/',length(starts),'windows\n')
}},finally={if(pdf_device%in%dev.list())dev.off(pdf_device)})
self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
checker<-file.path(dirname(normalizePath(self)),'09_pdf_qc.py')
python<-Sys.getenv('CRU_PYTHON',Sys.which('python3'))
if(!nzchar(python))stop('Python is required for PDF file QC')
info<-system2(python,c('-B',shQuote(checker),'--pdf',shQuote(temporary),
                      '--expected',as.character(length(starts))),stdout=TRUE)
status<-attr(info,'status')
if(!is.null(status)&&status!=0)stop('PDF verification failed')
page_qc<-read.delim(text=paste(info,collapse='\n'),na.strings='NA',check.names=FALSE)
pages<-page_qc$actual_pages
stopifnot(nrow(page_qc)==1,page_qc$page_count_status%in%c('verified','unavailable'),
          is.na(pages)||pages==length(starts),identical(tools::md5sum(paths),before))
stopifnot(file.rename(temporary,atlas))
index<-do.call(rbind,index);stopifnot(all(file.exists(index$png)))
write.table(index,file.path(out,'01_window_index.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(chrom=chrom,pdf_pages=pages,page_count_status=page_qc$page_count_status,
  page_count_backend=page_qc$page_count_backend,expected_windows=length(starts),png_n=nrow(index),
  coverage_start=min(index$start),coverage_end=max(index$end),window_bp=window,font=font),
  file.path(out,'02_render_QC.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(path=paths,md5=unname(before)),file.path(out,'03_plot_input_MD5.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PDF:',atlas,'\nPages:',pages,'\nPNG:',png_dir,'\n')
