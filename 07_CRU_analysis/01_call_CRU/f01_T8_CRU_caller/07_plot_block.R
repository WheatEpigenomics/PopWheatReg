#!/usr/bin/env Rscript
## fixed: character factor levels so the 100/200-kb strata are labelled, not NA
##*## Observed/random distance-stratified block enrichment, R vector output
a<-commandArgs(trailingOnly=TRUE);arg<-function(key,default=NULL){i<-match(key,a);if(is.na(i))return(default);a[i+1L]}
root<-arg('--analysis-dir');font<-arg('--font-family','Liberation Sans');line_scale<-as.numeric(arg('--line-scale','1.3'))
stopifnot(!is.null(root),file.exists(file.path(root,'07_ENRICHMENT.DONE.json')))
suppressPackageStartupMessages(library(ggplot2))
if(!nzchar(Sys.which('fc-match')))stop('fc-match required to verify font')
matched<-system2('fc-match',c('-f',shQuote('%{family}'),shQuote(font)),stdout=TRUE)
if(!any(grepl(tolower(font),tolower(matched),fixed=TRUE)))stop('Requested font not installed: ',font)
out<-file.path(root,'f01_figures');dir.create(out,showWarnings=FALSE)
input<-file.path(root,'02_enrichment_summary.tsv');data<-read.delim(input,check.names=FALSE)
d<-data[data$scope=='All_selected',];stopifnot(nrow(d)>0)
signature<-list(input_md5=unname(tools::md5sum(input)),font=font,line_scale=line_scale,
                script_md5=unname(tools::md5sum(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))))
stamp<-file.path(out,'03_PLOT.DONE.rds')
if(file.exists(stamp)) {
  if(identical(readRDS(stamp),signature)&&all(file.exists(file.path(out,c('01_block_boundary_enrichment.pdf','01_block_boundary_enrichment.png','02_all_final_enrichment.pdf','02_all_final_enrichment.png')))))quit(status=0)
  stop('Existing block figure options differ; preserve them in a new directory')
}
levels<-intersect(c('outer','internal','all_final'),unique(d$boundary_group));colours<-c(outer='#96CEDD',internal='#E89B6C',all_final='#B086CF')
labels<-setNames(vapply(levels,function(g){n<-unique(d$boundary_n[d$boundary_group==g]);paste(switch(g,outer='Outer boundaries',internal='Internal cuts',all_final='All-final boundaries'),paste0('(n=',n,')'))},character(1)),levels)
d$boundary_group<-factor(d$boundary_group,levels=levels)
observed<-data.frame(d,series='Observed',value=d$observed_fraction*100)
random<-data.frame(d,series='Random',value=d$random_fraction_mean*100)
plotdata<-rbind(observed,random);plotdata$window<-factor(as.character(plotdata$cutoff_bp),levels=c('25000','50000','100000','200000'),labels=c('\u226425','\u226450','\u2264100','\u2264200'))
d$window<-factor(as.character(d$cutoff_bp),levels=c('25000','50000','100000','200000'),labels=c('\u226425','\u226450','\u2264100','\u2264200'))
theme<-theme_classic(base_size=12,base_family=font)+theme(panel.grid=element_blank(),text=element_text(colour='#222222'),
  axis.line=element_line(linewidth=.7*line_scale),axis.ticks=element_line(linewidth=.7*line_scale),
  legend.position='top',legend.title=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))
make<-function(groups) {
  dd<-plotdata[plotdata$boundary_group%in%groups,];ann<-d[d$boundary_group%in%groups,]
  padding<-max(dd$value)*.07
  ggplot(dd,aes(window,value,group=series,colour=boundary_group,linetype=series,shape=series))+
    geom_errorbar(data=dd[dd$series=='Random',],aes(ymin=random_fraction_p025*100,ymax=random_fraction_p975*100),width=.08,linewidth=.45*line_scale,alpha=.45)+
    geom_line(linewidth=.75*line_scale)+geom_point(size=2.7)+
    geom_text(data=ann,aes(x=window,y=observed_fraction*100+padding,label=sprintf('%.2f\u00d7',observed_over_random)),
              inherit.aes=FALSE,colour='#222222',family=font,size=3.6)+
    scale_colour_manual(values=colours,guide='none')+scale_linetype_manual(values=c(Observed='solid',Random='dashed'))+
    scale_shape_manual(values=c(Observed=16,Random=0))+
    facet_wrap(~boundary_group,nrow=1,labeller=as_labeller(labels))+
    scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.15)))+
    labs(x='Distance to haplotype-block boundary (kb)',y='Boundary proportion (%)')+theme
}
savepair<-function(stem,p,width) {
  pdf<-file.path(out,paste0(stem,'.pdf'));pngfile<-file.path(out,paste0(stem,'.png'))
  if(file.exists(pdf)||file.exists(pngfile))stop('Figure exists without completed checkpoint')
  cairo_pdf(pdf,width=width,height=4.4,family=font);print(p);dev.off()
  png(pngfile,width=width,height=4.4,units='in',res=150,type='cairo',family=font);print(p);dev.off()
}
savepair('01_block_boundary_enrichment',make(levels),12)
savepair('02_all_final_enrichment',make('all_final'),5.3)
saveRDS(signature,stamp)
