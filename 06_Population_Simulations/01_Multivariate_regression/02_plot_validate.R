#!/usr/bin/env Rscript
# Draw the six original figures and independently validate all model results.
script_file <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
# Rscript may encode spaces as ~+~ in --file (including R 3.6).
if(!file.exists(script_file)) script_file <- gsub("~+~"," ",script_file,fixed=TRUE)
source(file.path(dirname(normalizePath(script_file)),"01_fit_models.R"))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L || !args[1] %in% c("all","plot","validate"))
  stop("Usage: Rscript 02_plot_validate.R {all|plot|validate} OUTPUT_DIR")
root <- normalizePath(args[2],mustWork=TRUE)


# plot_primary
plot_primary <- function(root) {
  figdir <- file.path(root,"figures")
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  source_data <- read_tsv(file.path(root,"plot_data/figure_source_data.tsv"))
  panel <- function(figure,panel_name) {
    d <- source_data[source_data$figure==figure & source_data$panel==panel_name, ]
    if(!nrow(d)) stop("Missing source data: ",figure,"/",panel_name)
    d[order(d$display_order), ]
  }
  COL <- c(agLink="#CC6677",genomic="#4477AA",TE="#228833",other="#999999")
  TEXT <- "#202020"; FONT_PT <- 7; LINE_PT <- .75
  LINE_LWD <- LINE_PT/(72/96) # Base-R lwd=1 is exactly 0.75 PostScript pt.
  fmt <- function(v) if(abs(v)<1e-4 && v!=0) format(v,scientific=TRUE,digits=2) else sprintf("%.4f",v)
  fmtp <- function(p) if(!is.finite(p)) "NA" else if(p<1e-4) format(p,scientific=TRUE,digits=2) else sprintf("%.4f",p)
  theme <- function() par(lwd=LINE_LWD,las=1,fg=TEXT,col=TEXT,family="Arial",cex=1)
  ax <- function(side,...) axis(side,...,lwd=LINE_LWD,lwd.ticks=LINE_LWD,tck=-.025,
                                 cex.axis=1,col=TEXT,col.axis=TEXT)
  ptitle <- function(title) mtext(title,side=3,adj=0,line=.25,font=2,cex=1)
  pdf_plot <- function(filename,width,height,draw) {
    cairo_pdf(file.path(figdir,filename),width=width,height=height,pointsize=FONT_PT,family="Arial")
    on.exit(dev.off())
    theme(); draw()
  }
  horizontal <- function(d,title,grouped=FALSE,mar=c(5,12.5,2,.8),footnote=NULL) {
    par(mar=mar,mgp=c(2.35,.65,0))
    d <- d[order(d$display_order,decreasing=TRUE), ]
    labels <- sub(" (","\n(",d$component,fixed=TRUE)
    xx <- d$unique_delta_R2
    xmax <- max(xx)*1.34; if(!is.finite(xmax)||xmax<=0) xmax <- .01
    pos <- barplot(xx,names.arg=labels,horiz=TRUE,col=unname(COL[d$colour_class]),
      border=TEXT,lwd=LINE_LWD,xlim=c(0,xmax),axes=FALSE,cex.names=1)
    ax(1); box(bty="l",lwd=LINE_LWD)
    title(xlab=if(grouped) expression("Grouped unique "*Delta*R^2) else expression("Unique "*Delta*R^2),
          line=2.4,cex.lab=1)
    text(xx+xmax*.018,pos,vapply(xx,fmt,character(1)),pos=4,cex=1)
    if(!is.null(footnote)) mtext(footnote,side=1,line=4,cex=1)
    ptitle(title)
  }

  pdf_plot("Fig_ABC_multivariate_confounder_control.pdf",9,3.4,function() {
    layout(matrix(1:3,nrow=1),widths=c(1.02,1.30,1.52)); theme()
    A <- panel("ABC","A")
    par(mar=c(5.2,4.6,2,.8),mgp=c(2.45,.65,0))
    vals <- A$R2; ymax <- max(vals)*1.42
    bp <- barplot(vals,names.arg=rep("",2),col=unname(COL[A$colour_class]),
      border=TEXT,lwd=LINE_LWD,ylim=c(0,ymax),axes=FALSE)
    mtext(sub(" ","\n",A$component,fixed=TRUE),side=1,at=bp,line=1.3,cex=1)
    ax(2); box(bty="l",lwd=LINE_LWD)
    title(ylab=expression("Model "*R^2),line=2.55,cex.lab=1)
    text(bp,vals+ymax*.025,sprintf("%.3f",vals),cex=1)
    yb <- ymax*.79
    segments(bp[1],yb,bp[2],yb,lwd=LINE_LWD)
    segments(bp,yb,bp,yb-ymax*.025,lwd=LINE_LWD)
    text(mean(bp),yb+ymax*.105,paste0("+agLinks: ΔR²=",fmt(A$step_delta_R2[2]),
      "\nLRT P=",fmtp(A$LRT_P[2]),"; wild P=",fmtp(A$wild_bootstrap_P[2])),cex=1)
    mtext(paste0("n=",A$n_blocks[1]," blocks"),side=1,line=4.25,cex=1)
    ptitle("A  Model improvement")
    horizontal(panel("ABC","B"),"B  Drop-one contribution",mar=c(5,9.2,2,.8),
               footnote="Predictor removed from the same full model")
    horizontal(panel("ABC","C"),"C  Regulatory vs genomic contribution",grouped=TRUE,
               footnote="Grouped drop-one; chromosome retained")
  })

  pdf_plot("Fig_4Mb_context_ablation.pdf",4.5,3.1,function() {
    d <- panel("4Mb","ablation")
    horizontal(d,"4-Mb context ablation",mar=c(5,12.2,2,1.1),
      footnote=paste0("Drop-one from full model; n=",d$n_blocks[1]," blocks"))
  })
  pdf_plot("Fig_E_regulatory_vs_genomic_contribution.pdf",4.5,3.55,function() {
    horizontal(panel("E","E"),"E  Regulatory vs genomic contribution",grouped=TRUE,
      mar=c(5,13.2,2,.9),footnote="Grouped drop-one; chromosome retained")
  })
  pdf_plot("Fig_F_sequential_model_improvement.pdf",4.8,4.1,function() {
    d <- panel("F","F")
    par(mar=c(8,4.6,2.2,1),mgp=c(2.45,.65,0))
    ymax <- max(d$R2)*1.72
    bp <- barplot(d$R2,names.arg=rep("",3),col=unname(COL[d$colour_class]),
      border=TEXT,lwd=LINE_LWD,ylim=c(0,ymax),axes=FALSE)
    ax(2); box(bty="l",lwd=LINE_LWD)
    title(ylab=expression("Model "*R^2),line=2.55,cex.lab=1)
    mtext(c("baseline","baseline + TE","baseline + TE"),side=1,at=bp,line=1.2,cex=1)
    mtext("+ agLink",side=1,at=bp[3],line=2.4,cex=1)
    text(bp,d$R2+ymax*.025,sprintf("%.4f",d$R2),cex=1)
    for(k in 2:3) {
      yb <- ymax*c(.70,.88)[k-1]
      segments(bp[k-1],yb,bp[k],yb,lwd=LINE_LWD)
      segments(bp[c(k-1,k)],yb,bp[c(k-1,k)],yb-ymax*.025,lwd=LINE_LWD)
      text(mean(bp[c(k-1,k)]),yb+ymax*.063,
        paste0(if(k==2) "+TE" else "+agLinks",": ΔR²=",fmt(d$step_delta_R2[k]),
          "\nLRT P=",fmtp(d$LRT_P[k]),"; wild P=",fmtp(d$wild_bootstrap_P[k])),cex=1)
    }
    ptitle("F  Sequential model improvement")
    mtext(paste0("n=",d$n_blocks[1]," blocks; SNP density excluded"),side=1,line=3.9,cex=1)
    mtext("Baseline: Gene density + Recombination",side=1,line=5.3,cex=1)
    mtext("+ ACR density + chromosome",side=1,line=6.6,cex=1)
  })
  cat("DONE: four PDF figures generated from the single source-data table.\n")
}


# plot_snp_sequential
plot_snp_sequential <- function(root) {
  figdir <- file.path(root,"figures")
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  d <- read_tsv(file.path(root,"plot_data/Fig_F_sequential_with_SNP.tsv"))
  d <- d[order(d$display_order), ]
  stopifnot(nrow(d)==3L,all(d$SNP_density_used),all(is.finite(d$R2)))
  COL <- c(agLink="#CC6677",genomic="#4477AA",TE="#228833",other="#999999")
  TEXT <- "#202020"; FONT_PT <- 7; LINE_PT <- .75; LWD <- LINE_PT/(72/96)
  fmt <- function(v) if(abs(v)<1e-4 && v!=0) format(v,scientific=TRUE,digits=2) else sprintf("%.4f",v)
  fmtp <- function(p) if(!is.finite(p)) "NA" else if(p<1e-4) format(p,scientific=TRUE,digits=2) else sprintf("%.4f",p)
  cairo_pdf(file.path(figdir,"Fig_F_sequential_model_improvement_with_SNP.pdf"),
            width=4.8,height=4.1,pointsize=FONT_PT,family="Arial")
  par(lwd=LWD,las=1,fg=TEXT,col=TEXT,family="Arial",cex=1,
      mar=c(8,4.6,2.2,1),mgp=c(2.45,.65,0))
  ymax <- max(d$R2)*1.72
  bp <- barplot(d$R2,names.arg=rep("",3),col=unname(COL[d$colour_class]),
                border=TEXT,lwd=LWD,ylim=c(0,ymax),axes=FALSE)
  axis(2,lwd=LWD,lwd.ticks=LWD,tck=-.025,cex.axis=1,col=TEXT,col.axis=TEXT)
  box(bty="l",lwd=LWD)
  title(ylab=expression("Model "*R^2),line=2.55,cex.lab=1)
  mtext(c("baseline","baseline + TE","baseline + TE"),side=1,at=bp,line=1.2,cex=1)
  mtext("+ agLink",side=1,at=bp[3],line=2.4,cex=1)
  text(bp,d$R2+ymax*.025,sprintf("%.4f",d$R2),cex=1)
  for(k in 2:3) {
    yb <- ymax*c(.70,.88)[k-1]
    segments(bp[k-1],yb,bp[k],yb,lwd=LWD)
    segments(bp[c(k-1,k)],yb,bp[c(k-1,k)],yb-ymax*.025,lwd=LWD)
    text(mean(bp[c(k-1,k)]),yb+ymax*.063,
      paste0(if(k==2) "+TE" else "+agLinks",": ΔR²=",fmt(d$step_delta_R2[k]),
        "\nLRT P=",fmtp(d$LRT_P[k]),"; wild P=",fmtp(d$wild_bootstrap_P[k])),cex=1)
  }
  mtext("F  Sequential model improvement (with SNP)",side=3,adj=0,line=.25,font=2,cex=1)
  mtext(paste0("n=",d$n_blocks[1]," blocks; SNP density included"),side=1,line=3.9,cex=1)
  mtext("Baseline: Gene density + Recombination",side=1,line=5.3,cex=1)
  mtext("+ ACR density + SNP density + chromosome",side=1,line=6.6,cex=1)
  invisible(dev.off())
  cat("DONE: additional F PDF with SNP density, 7-pt Arial and 0.75-pt lines.\n")
}


# plot_snp_drop_one
plot_snp_drop_one <- function(root) {
  figdir <- file.path(root,"figures")
  d <- read_tsv(file.path(root,"plot_data/Fig_B_drop_one_with_SNP.tsv"))
  d <- d[order(d$display_order), ]
  stopifnot(nrow(d)==6L,identical(d$display_order,1:6),all(d$SNP_density_in_full_model),
            all(is.finite(d$unique_delta_R2)),all(d$unique_delta_R2>=0))
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  FONT_PT <- 7; LINE_PT <- .75; LWD <- LINE_PT/(72/96)
  TEXT <- "#202020"
  cairo_pdf(file.path(figdir,"Fig_B_drop_one_contribution_with_SNP.pdf"),
            width=2.65,height=2.9,pointsize=FONT_PT,family="Arial")
  par(mai=c(.42,1.02,.50,.12),mgp=c(1.9,.45,0),family="Arial",
      cex=1,cex.axis=1,cex.lab=1,las=1,lwd=LWD,fg=TEXT,col=TEXT,xaxs="i",yaxs="i")
  xmax <- max(.065,max(d$unique_delta_R2)*1.32)
  y <- rev(seq_len(nrow(d)))
  plot.new(); plot.window(xlim=c(0,xmax),ylim=c(.3,6.7))
  rect(0,y-.40,d$unique_delta_R2,y+.40,col=d$colour_hex,border=TEXT,lwd=LWD)
  labels <- ifelse(d$component=="Recombination intensity","Recombination\nintensity",d$component)
  axis(2,at=y,labels=labels,tick=FALSE,line=0,hadj=1,cex.axis=1)
  axis(1,at=c(0,.03,.06),labels=c("0","0.03","0.06"),tck=-.025,
       lwd=LWD,lwd.ticks=LWD,cex.axis=1)
  box(bty="l",lwd=LWD)
  value_labels <- sprintf("%.4f",d$unique_delta_R2)
  label_x <- d$unique_delta_R2+xmax*.022
  # Fail rather than silently clip values when data or fonts change.
  stopifnot(all(label_x+strwidth(value_labels,units="user")<xmax))
  text(label_x,y,value_labels,adj=c(0,.5),cex=1)
  mtext("Unique ΔR²",side=1,line=1.65,cex=1)
  mtext("Predictor removed from\nthe same full model",side=3,line=.65,cex=1)
  invisible(dev.off())
  cat("DONE: six-bar SNP-included drop-one PDF, Arial 7 pt, lines 0.75 pt.\n")
}


# validate_primary
validate_primary <- function(root) {
  expected_n <- expected_blocks(root)
  report <- list()
  check <- function(ok,label,detail="") {
    report[[length(report)+1L]] <<- data.frame(check=label,status=if(isTRUE(ok)) "PASS" else "FAIL",
                                             detail=detail,stringsAsFactors=FALSE)
    if(!isTRUE(ok)) {
      write_tsv(do.call(rbind,report),file.path(root,"docs/validation_report.tsv"))
      stop("Validation failed: ",label," ",detail)
    }
  }
  near <- function(a,b,label,tol=1e-9) {
    check(length(a)==length(b) && all(is.finite(a)) && all(is.finite(b)) &&
          max(abs(a-b),0)<=tol,label,paste("maximum allowed absolute error",tol))
  }
  check_p <- function(a,b,label) {
    same_zero <- (a==0 & b==0)
    ok <- is.finite(a) & is.finite(b) & a>=0 & b>=0 & a<=1 & b<=1
    check(all(ok) && all(same_zero | (a>0 & b>0 & abs(log(pmax(a,.Machine$double.xmin))-
          log(pmax(b,.Machine$double.xmin)))<1e-7)),label)
  }
  r2 <- function(m) summary(m)$r.squared
  partial <- function(m0,m1) 1-sum(resid(m1)^2)/sum(resid(m0)^2)
  lrt <- function(m0,m1) {
    chi <- 2*(as.numeric(logLik(m1))-as.numeric(logLik(m0)))
    df <- attr(logLik(m1),"df")-attr(logLik(m0),"df")
    c(chi=chi,df=df,p=pchisq(chi,df,lower.tail=FALSE))
  }
  # Independent cluster-score implementation, using rowsum instead of fit_models.R's loop.
  cluster_vcov <- function(m,cluster) {
    X <- model.matrix(m); n <- nrow(X); p <- ncol(X); G <- nlevels(factor(cluster))
    score <- rowsum(X*as.numeric(resid(m)),cluster,reorder=TRUE)
    bread <- solve(crossprod(X))
    (G/(G-1))*((n-1)/(n-p))*bread %*% crossprod(score) %*% bread
  }
  cluster_p <- function(m,cluster,terms) {
    V <- cluster_vcov(m,cluster); b <- coef(m)[terms]; G <- nlevels(factor(cluster))
    value <- if(length(terms)==1L) 2*pt(abs(b/sqrt(V[terms,terms])),G-1,lower.tail=FALSE) else
      pf(as.numeric(crossprod(b,solve(V[terms,terms,drop=FALSE],b)))/length(terms),
         length(terms),G-1,lower.tail=FALSE)
    unname(value)
  }

  input <- read_tsv(file.path(root,"data/analysis_blocks.tsv.gz"))
  nested <- read_tsv(file.path(root,"results/nested_models.tsv"))
  coefficients <- read_tsv(file.path(root,"results/model_coefficients.tsv"))
  contributions <- read_tsv(file.path(root,"results/contributions.tsv"))
  sensitivity <- read_tsv(file.path(root,"results/sensitivity.tsv"))
  diagnostics <- read_tsv(file.path(root,"results/diagnostics.tsv"))
  nulls <- read_tsv(file.path(root,"results/bootstrap_null.tsv.gz"))
  source_data <- read_tsv(file.path(root,"plot_data/figure_source_data.tsv"))
  settings <- read_tsv(file.path(root,"docs/analysis_settings.tsv"))
  setting <- function(key) settings$value[match(key,settings$setting)]
  check(identical(names(input),c("context",input_columns)),"input_schema")
  check(!any(grepl("snp",c(names(input),coefficients$term,nested$formula,contributions$component),ignore.case=TRUE)),
        "SNP_excluded_from_input_and_models")
  check(!anyDuplicated(paste(input$context,input$block_id)),"unique_context_block_identifiers")
  check(setequal(unique(input$context),context_order),"two_contexts_present")
  check(as.integer(setting("bootstrap_replicates"))==9999L,"9999_bootstrap_replicates")
  check(!anyDuplicated(settings$setting) && is.finite(as.integer(setting("random_seed"))) &&
        identical(as.character(c(setting("RNG_uniform"),setting("RNG_normal"),setting("RNG_sample"))),
                  c("Mersenne-Twister","Inversion","Rejection")),"recorded_RNG_settings")
  check(nrow(nested)==6L && nrow(coefficients)==10L && nrow(contributions)==24L &&
        nrow(sensitivity)==8L && nrow(nulls)==39996L,"consolidated_table_row_counts")
  check(identical(input$block_id[input$context=="2Mb"],input$block_id[input$context=="4Mb"]),
        "identical_ordered_blocks_between_contexts")
  for(ctx in context_order) {
    prefix <- paste0(ctx,": ")
    x <- input[input$context==ctx, ]; x$chrom <- factor(x$chrom)
    check(nrow(x)==as.integer(setting("frozen_blocks_per_context")),paste0(prefix,"frozen_sample_size"))
    check(all(is.finite(as.matrix(x[,setdiff(input_columns,c("block_id","block_name","chrom"))]))),
          paste0(prefix,"finite_input_values"))
    check(all(x$length_bp==x$end-x$start & x$length_bp>0) && all(x$start>=0) &&
          all(x$context_length_bp==x$context_end-x$context_start & x$context_length_bp>0) &&
          all(x$context_start>=0 & x$context_length_bp<=if(ctx=="2Mb") 2000000 else 4000000),
          paste0(prefix,"block_and_context_coordinates"))
    check(all(x$tested_gene_count>0 & x$aglinks_per_tested_gene>=0 & x$aglinks_same_block_per_tested_gene>=0 &
          x$context_gene_density_per_mb>=0 & x$context_acr_density_per_mb>=0 &
          x$context_recombination_rate>=0 & x$context_recomb_coverage_fraction>=.5 &
          x$context_recomb_coverage_fraction<=1 & x$context_te_content>=0 & x$context_te_content<=1),
          paste0(prefix,"input_ranges_and_coverage"))
    near(x$aglinks_per_tested_gene,x$aglink_count_gene_block/x$tested_gene_count,paste0(prefix,"agLink_ratio"))
    z <- function(v) as.numeric(scale(v))
    x$y <- log1p(x$length_bp); x$ag_z <- z(log1p(x$aglinks_per_tested_gene))
    x$gene_z <- z(log1p(x$context_gene_density_per_mb)); x$recomb_z <- z(x$context_recombination_rate)
    x$te_z <- z(x$context_te_content); x$acr_z <- z(log1p(x$context_acr_density_per_mb))
    models <- list(lm(y~gene_z+recomb_z+acr_z+chrom,data=x),
                   lm(y~gene_z+recomb_z+te_z+acr_z+chrom,data=x),
                   lm(y~ag_z+gene_z+recomb_z+te_z+acr_z+chrom,data=x))
    full <- models[[3]]; reduced <- models[[2]]
    check(full$rank==ncol(model.matrix(full)),paste0(prefix,"full_model_rank"))
    tab <- nested[nested$context==ctx, ]; tab <- tab[order(tab$stage), ]
    near(tab$R2,vapply(models,r2,numeric(1)),paste0(prefix,"nested_model_R2"))
    near(tab$adjusted_R2,vapply(models,function(m) summary(m)$adj.r.squared,numeric(1)),
         paste0(prefix,"nested_model_adjusted_R2"))
    check(all(tab$n_blocks==nrow(x)) && all(is.na(tab$step_delta_R2[1])) &&
          all(is.na(tab$LRT_P[1])) && all(is.na(tab$wild_bootstrap_P[1])),paste0(prefix,"nested_model_baseline"))
    near(tab$step_delta_R2[2:3],diff(tab$R2),paste0(prefix,"sequential_delta_R2"))
    for(k in 2:3) {
      ll <- lrt(models[[k-1]],models[[k]])
      near(tab$LRT_chisq[k],ll["chi"],paste0(prefix,"step_",k,"_LRT_chisq"))
      check(tab$LRT_df[k]==ll["df"],paste0(prefix,"step_",k,"_LRT_df"))
      check_p(tab$LRT_P[k],ll["p"],paste0(prefix,"step_",k,"_LRT_P"))
      near(tab$step_partial_R2[k],partial(models[[k-1]],models[[k]]),paste0(prefix,"step_",k,"_partial_R2"))
      term <- if(k==2) "te_z" else "ag_z"
      near(tab$chromosome_cluster_beta[k],unname(coef(models[[k]])[term]),paste0(prefix,"step_",k,"_beta"))
      check_p(tab$chromosome_cluster_P[k],cluster_p(models[[k]],x$chrom,term),paste0(prefix,"step_",k,"_cluster_P"))
      nn <- nulls[nulls$context==ctx & nulls$comparison==if(k==2) "add_TE" else "add_agLink", ]
      check(nrow(nn)==9999L && identical(nn$iteration,seq_len(9999L)) && all(is.finite(nn$null_cluster_t)),
            paste0(prefix,"step_",k,"_bootstrap_null"))
      tt <- unname(coef(models[[k]])[term])/sqrt(cluster_vcov(models[[k]],x$chrom)[term,term])
      near(tab$wild_bootstrap_P[k],(1+sum(abs(nn$null_cluster_t)>=abs(tt)))/10000,
           paste0(prefix,"step_",k,"_bootstrap_P"),tol=1e-12)
    }
    cc <- coefficients[coefficients$context==ctx, ]
    check(setequal(cc$term,predictor_terms),paste0(prefix,"coefficient_terms"))
    near(cc$beta,unname(coef(full)[cc$term]),paste0(prefix,"coefficients"))
    ols <- coef(summary(full)); V <- cluster_vcov(full,x$chrom); G <- nlevels(x$chrom)
    near(cc$OLS_SE,ols[cc$term,"Std. Error"],paste0(prefix,"OLS_SE"))
    check_p(cc$OLS_P,ols[cc$term,"Pr(>|t|)"],paste0(prefix,"OLS_P"))
    se <- sqrt(diag(V)[cc$term])
    near(cc$chromosome_cluster_SE,se,paste0(prefix,"cluster_SE"))
    check_p(cc$chromosome_cluster_P,vapply(cc$term,function(t) cluster_p(full,x$chrom,t),numeric(1)),
            paste0(prefix,"cluster_coefficient_P"))
    near(cc$cluster_CI_low,cc$beta-qt(.975,G-1)*se,paste0(prefix,"cluster_95CI_low"))
    near(cc$cluster_CI_high,cc$beta+qt(.975,G-1)*se,paste0(prefix,"cluster_95CI_high"))
    for(decomp in names(contribution_groups)) {
      ct <- contributions[contributions$context==ctx & contributions$decomposition==decomp, ]
      groups <- contribution_groups[[decomp]]
      check(setequal(ct$component,names(groups)) && all(ct$n_blocks==nrow(x)),paste0(prefix,decomp,"_groups"))
      for(i in seq_len(nrow(ct))) {
        mr <- lm(reformulate(c(setdiff(predictor_terms,groups[[ct$component[i]]]),"chrom"),response="y"),data=x)
        label <- paste0(prefix,decomp,"_",i)
        near(c(ct$full_R2[i],ct$reduced_R2[i],ct$unique_delta_R2[i],ct$partial_R2[i]),
             c(r2(full),r2(mr),r2(full)-r2(mr),partial(mr,full)),paste0(label,"_R2"))
        ll <- lrt(mr,full)
        near(ct$LRT_chisq[i],ll["chi"],paste0(label,"_LRT_chisq"))
        check(ct$LRT_df[i]==ll["df"],paste0(label,"_LRT_df"))
        check_p(ct$LRT_P[i],ll["p"],paste0(label,"_LRT_P"))
      }
    }
    x$has_aglink <- as.integer(x$aglinks_per_tested_gene>0); pos <- x$has_aglink==1
    x$positive_log_degree_centered <- 0
    degree_center <- median(log(x$aglinks_per_tested_gene[pos]))
    x$positive_log_degree_centered[pos] <- log(x$aglinks_per_tested_gene[pos])-degree_center
    two <- lm(y~has_aglink+positive_log_degree_centered+gene_z+recomb_z+te_z+acr_z+chrom,data=x)
    x$ag_same_z <- z(log1p(x$aglinks_same_block_per_tested_gene))
    same <- lm(y~ag_same_z+gene_z+recomb_z+te_z+acr_z+chrom,data=x)
    xs <- x[x$candidate_acr_gene_pair_proxy>0 & x$tested_gene_count>0, ]
    xs$log_tested <- log(xs$tested_gene_count); xs$log_pairs <- log(xs$candidate_acr_gene_pair_proxy)
    pp <- xs$has_aglink==1; xs$positive_log_degree_centered <- 0
    xs$positive_log_degree_centered[pp] <- log(xs$aglinks_per_tested_gene[pp])-median(log(xs$aglinks_per_tested_gene[pp]))
    strict0 <- lm(y~splines::ns(log_tested,df=4)+splines::ns(log_pairs,df=3)+gene_z+recomb_z+te_z+acr_z+chrom,data=xs)
    strict1 <- lm(y~has_aglink+positive_log_degree_centered+splines::ns(log_tested,df=4)+
                   splines::ns(log_pairs,df=3)+gene_z+recomb_z+te_z+acr_z+chrom,data=xs)
    st <- sensitivity[sensitivity$context==ctx, ]
    near(st$delta_R2,c(r2(full)-r2(reduced),r2(two)-r2(reduced),r2(same)-r2(reduced),r2(strict1)-r2(strict0)),
         paste0(prefix,"sensitivity_delta_R2"))
    near(st$partial_R2,c(partial(reduced,full),partial(reduced,two),partial(reduced,same),partial(strict0,strict1)),
         paste0(prefix,"sensitivity_partial_R2"))
    check_p(st$chromosome_cluster_P,c(cluster_p(full,x$chrom,"ag_z"),
      cluster_p(two,x$chrom,c("has_aglink","positive_log_degree_centered")),cluster_p(same,x$chrom,"ag_same_z"),
      cluster_p(strict1,xs$chrom,c("has_aglink","positive_log_degree_centered"))),paste0(prefix,"sensitivity_cluster_P"))
    check(all(st$wild_cluster_P>=1e-4 & st$wild_cluster_P<=1),paste0(prefix,"sensitivity_bootstrap_P_range"))
    diag <- diagnostics[diagnostics$context==ctx, ]
    vf <- diag[diag$diagnostic=="VIF", ]
    expected_vif <- vapply(vf$variable_1,function(t) {
      other <- c(setdiff(predictor_terms,t),"chrom")
      mm <- lm(reformulate(other,response=t),data=x)
      1/(1-r2(mm))
    },numeric(1))
    near(vf$value,expected_vif,paste0(prefix,"VIF"))
    cor_data <- data.frame(block_length=x$length_bp,aglinks_per_tested_gene=x$aglinks_per_tested_gene,
      gene_density=x$context_gene_density_per_mb,recombination=x$context_recombination_rate,TE=x$context_te_content,
      ACR_density=x$context_acr_density_per_mb,tested_gene_count=x$tested_gene_count,candidate_pair_proxy=x$candidate_acr_gene_pair_proxy)
    cm <- cor(cor_data,method="spearman"); sp <- diag[diag$diagnostic=="Spearman", ]
    near(sp$value,cm[cbind(match(sp$variable_1,rownames(cm)),match(sp$variable_2,colnames(cm)))],paste0(prefix,"Spearman"))
    qc <- diag[diag$diagnostic=="analysis_QC", ]
    expected <- c(n_blocks=nrow(x),n_chromosomes=G,primary_reduced_R2=r2(reduced),primary_full_R2=r2(full),
      primary_agLink_delta_R2=r2(full)-r2(reduced),primary_agLink_partial_R2=partial(reduced,full),
      primary_LRT_P=lrt(reduced,full)[["p"]],primary_chr_cluster_P=cluster_p(full,x$chrom,"ag_z"),
      primary_wild_cluster_P=tab$wild_bootstrap_P[3],two_part_degree_center=exp(degree_center))
    near(qc$value,unname(expected[qc$variable_1]),paste0(prefix,"analysis_QC"))
  }
  expected_source <- make_figure_source(nested,contributions)
  check(isTRUE(all.equal(source_data,expected_source,tolerance=1e-12,check.attributes=FALSE)),
        "every_figure_value_matches_statistical_results")
  check(nrow(source_data)==22L && !anyDuplicated(paste(source_data$figure,source_data$panel,source_data$display_order)),
        "unique_figure_source_rows")
  for(f in c("01_fit_models.R","02_plot_validate.R")) {
    invisible(parse(file.path(dirname(normalizePath(script_file)),f)))
  }
  check(TRUE,"all_R_scripts_parse")
  dict <- read_tsv(file.path(root,"docs/data_dictionary.tsv"))
  check(all(names(input) %in% dict$field),"all_input_fields_documented")
  expected_pdfs <- c("Fig_ABC_multivariate_confounder_control.pdf","Fig_4Mb_context_ablation.pdf",
                     "Fig_E_regulatory_vs_genomic_contribution.pdf","Fig_F_sequential_model_improvement.pdf",
                     "Fig_F_sequential_model_improvement_with_SNP.pdf",
                     "Fig_B_drop_one_contribution_with_SNP.pdf")
  actual_files <- list.files(file.path(root,"figures"))
  check(setequal(actual_files,expected_pdfs),"exactly_six_PDF_figures_no_raster_or_temporary_files")
  for(f in expected_pdfs) {
    con <- file(file.path(root,"figures",f),"rb")
    header <- rawToChar(readBin(con,"raw",n=5)); close(con)
    check(header=="%PDF-" && file.info(file.path(root,"figures",f))$size>1000,
          paste0("valid_PDF_header: ",f))
  }
  write_tsv(do.call(rbind,report),file.path(root,"docs/validation_report.tsv"))
  cat("PASS:",length(report),"checks; frozen inputs, both models, sensitivities, bootstraps and all figure values agree.\n")
}


# validate_snp_sequential
validate_snp_sequential <- function(root) {
  expected_n <- expected_blocks(root)
  report <- list()
  check <- function(ok,label) {
    report[[length(report)+1L]] <<- data.frame(check=label,status=if(isTRUE(ok)) "PASS" else "FAIL")
    if(!isTRUE(ok)) {
      write_tsv(do.call(rbind,report),file.path(root,"docs/SNP_included_validation_report.tsv"))
      stop("SNP-included validation failed: ",label)
    }
  }
  near <- function(a,b,label,tol=1e-9) check(length(a)==length(b) && all(is.finite(a)) &&
    all(is.finite(b)) && max(abs(a-b),0)<tol,label)
  check_p <- function(a,b,label) check(length(a)==length(b) && all(is.finite(a)) && all(is.finite(b)) &&
    all(a>0 & a<=1 & b>0 & b<=1) && max(abs(log(a)-log(b)),0)<1e-7,label)
  input <- read_tsv(file.path(root,"data/analysis_blocks.tsv.gz"))
  x <- input[input$context=="2Mb", ]
  s <- read_tsv(file.path(root,"data/SNP_density_2Mb.tsv.gz"))
  tab <- read_tsv(file.path(root,"results/nested_models_with_SNP.tsv"))
  source_data <- read_tsv(file.path(root,"plot_data/Fig_F_sequential_with_SNP.tsv"))
  null <- read_tsv(file.path(root,"results/bootstrap_null_with_SNP.tsv.gz"))
  coeff <- read_tsv(file.path(root,"results/model_coefficients_with_SNP.tsv"))
  rows <- read_tsv(file.path(root,"results/model_rows_with_SNP.tsv.gz"))
  settings <- read_tsv(file.path(root,"docs/SNP_included_settings.tsv"))
  provenance <- read_tsv(file.path(root,"docs/SNP_covariate_provenance.tsv"))
  check(nrow(x)==expected_n && nrow(s)==expected_n && !anyDuplicated(s$block_id),"unique_SNP_rows_match_settings")
  check(identical(s$block_id,x$block_id),"same_original_ordered_block_IDs")
  for(col in c("context","chrom","start","end","context_start","context_end","context_length_bp"))
    check(identical(as.character(x[[col]]),as.character(s[[col]])),paste0("SNP_join_",col))
  check(all(is.finite(s$context_snp_density_per_mb)) && all(s$context_snp_density_per_mb>=0),"finite_nonnegative_SNP_density")
  check(all(s$context_snp_count>=0 & s$context_snp_count==round(s$context_snp_count)),"nonnegative_integer_SNP_counts")
  near(s$context_snp_density_per_mb,s$context_snp_count*1e6/s$context_length_bp,"SNP_density_from_count_and_context_length",1e-8)
  check(as.numeric(provenance$value[provenance$setting=="maximum_common_covariate_error"])==0,
        "original_archive_and_current_common_covariates_identical")
  x$context_snp_density_per_mb <- s$context_snp_density_per_mb
  x$context_snp_count <- s$context_snp_count
  x$chrom <- factor(x$chrom)
  x$y <- log1p(x$length_bp)
  x$ag_z <- as.numeric(scale(log1p(x$aglinks_per_tested_gene)))
  x$gene_z <- as.numeric(scale(log1p(x$context_gene_density_per_mb)))
  x$acr_z <- as.numeric(scale(log1p(x$context_acr_density_per_mb)))
  x$recomb_z <- as.numeric(scale(x$context_recombination_rate))
  x$te_z <- as.numeric(scale(x$context_te_content))
  x$snp_z <- as.numeric(scale(log1p(x$context_snp_density_per_mb)))
  check(nlevels(x$chrom)==21L,"21_original_chromosomes")
  check(identical(rows$block_id,x$block_id),"saved_model_rows_preserve_block_IDs")
  for(col in setdiff(names(rows),c("block_id","chrom"))) near(rows[[col]],x[[col]],paste0("model_input_",col))
  models <- list(lm(y~gene_z+recomb_z+snp_z+acr_z+chrom,data=x),
    lm(y~gene_z+recomb_z+te_z+snp_z+acr_z+chrom,data=x),
    lm(y~ag_z+gene_z+recomb_z+te_z+snp_z+acr_z+chrom,data=x))
  check(identical(tab$stage,1:3) && all(tab$n_blocks==expected_n) && all(tab$SNP_density_used),"three_SNP_included_stages")
  check(all(grepl("snp_z",tab$formula,fixed=TRUE)),"SNP_retained_in_ALL_three_models")
  check(all(vapply(models,function(m) m$rank==ncol(model.matrix(m)),logical(1))),"all_SNP_models_full_rank")
  near(tab$R2,vapply(models,function(m) summary(m)$r.squared,numeric(1)),"independent_SNP_model_R2")
  near(tab$adjusted_R2,vapply(models,function(m) summary(m)$adj.r.squared,numeric(1)),"independent_adjusted_R2")
  check(all(is.na(tab$step_delta_R2[1])) && is.na(tab$LRT_P[1]) && is.na(tab$wild_bootstrap_P[1]),"baseline_no_preceding_test")
  near(tab$step_delta_R2[2:3],diff(tab$R2),"adjacent_step_delta_R2")
  no_snp <- read_tsv(file.path(root,"results/nested_models.tsv"))
  no_snp <- no_snp[no_snp$context=="2Mb", ]; no_snp <- no_snp[order(no_snp$stage), ]
  check(all(tab$R2>=no_snp$R2-1e-12),"adding_SNP_does_not_reduce_in_sample_R2")
  cluster_vcov_independent <- function(X,u,cluster) {
    n <- nrow(X); p <- ncol(X); G <- nlevels(factor(cluster))
    scores <- rowsum(X*as.numeric(u),cluster,reorder=TRUE)
    bread <- solve(crossprod(X))
    (G/(G-1))*((n-1)/(n-p))*bread%*%crossprod(scores)%*%bread
  }
  for(k in 1:3) {
    m <- models[[k]]; X <- model.matrix(m); V <- cluster_vcov_independent(X,resid(m),x$chrom)
    cc <- coeff[coeff$stage==k, ]; ct <- coef(summary(m))
    terms <- rownames(ct)[rownames(ct) %in% c("ag_z","gene_z","recomb_z","te_z","snp_z","acr_z")]
    ct <- ct[terms,,drop=FALSE]
    check(identical(cc$term,terms),paste0("stage_",k,"_coefficient_terms"))
    near(cc$beta,unname(coef(m)[terms]),paste0("stage_",k,"_coefficients"))
    near(cc$OLS_SE,ct[,2],paste0("stage_",k,"_OLS_SE"))
    check_p(cc$OLS_P,ct[,4],paste0("stage_",k,"_OLS_P"))
    se <- sqrt(pmax(0,diag(V)[terms])); crit <- qt(.975,20)
    near(cc$chromosome_cluster_SE,unname(se),paste0("stage_",k,"_independent_cluster_SE"))
    near(cc$cluster_CI_low,unname(coef(m)[terms]-crit*se),paste0("stage_",k,"_cluster_CI_low"))
    near(cc$cluster_CI_high,unname(coef(m)[terms]+crit*se),paste0("stage_",k,"_cluster_CI_high"))
    if(k==1L) next
    term <- if(k==2) "te_z" else "ag_z"
    expected_seed <- if(k==2) 20260914L else 20260910L
    previous <- models[[k-1]]
    chi <- 2*(as.numeric(logLik(m))-as.numeric(logLik(previous)))
    df <- attr(logLik(m),"df")-attr(logLik(previous),"df")
    near(tab$LRT_chisq[k],chi,paste0("stage_",k,"_LRT_statistic"))
    check(tab$LRT_df[k]==df,paste0("stage_",k,"_LRT_df"))
    check_p(tab$LRT_P[k],pchisq(chi,df,lower.tail=FALSE),paste0("stage_",k,"_LRT_P"))
    near(tab$step_partial_R2[k],1-sum(resid(m)^2)/sum(resid(previous)^2),paste0("stage_",k,"_partial_R2"))
    tt <- unname(coef(m)[term]/sqrt(V[term,term]))
    near(tab$observed_cluster_t[k],tt,paste0("stage_",k,"_observed_cluster_t"))
    check_p(tab$chromosome_cluster_P[k],unname(2*pt(abs(tt),20,lower.tail=FALSE)),paste0("stage_",k,"_cluster_P"))
    nn <- null[null$comparison==if(k==2) "add_TE" else "add_agLink", ]
    check(nrow(nn)==9999L && identical(nn$iteration,seq_len(9999L)) && all(is.finite(nn$null_cluster_t)),
          paste0("stage_",k,"_9999_saved_null_statistics"))
    check(tab$wild_bootstrap_B[k]==9999L && tab$wild_bootstrap_seed[k]==expected_seed,paste0("stage_",k,"_original_B_and_seed"))
    near(tab$wild_bootstrap_P[k],(1+sum(abs(nn$null_cluster_t)>=abs(tt)))/10000,
         paste0("stage_",k,"_bootstrap_P_from_nulls"),1e-12)
    # Recreate selected draws independently with lm.fit and rowsum scores.
    RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(expected_seed)
    webb <- c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
    replay <- c(1L,2L,3L,10L,25L,100L,3333L,9999L); replayed <- numeric(length(replay))
    for(b in seq_len(9999L)) {
      w <- sample(webb,21,replace=TRUE)
      if(!b %in% replay) next
      ys <- fitted(previous)+resid(previous)*w[as.integer(x$chrom)]
      mm <- lm.fit(X,ys); VV <- cluster_vcov_independent(X,mm$residuals,x$chrom)
      replayed[match(b,replay)] <- mm$coefficients[term]/sqrt(VV[term,term])
    }
    near(replayed,nn$null_cluster_t[replay],paste0("stage_",k,"_independent_bootstrap_replay"),1e-8)
  }
  check(isTRUE(all.equal(source_data[,names(tab)],tab,tolerance=1e-12,check.attributes=FALSE)),"every_SNP_figure_value_matches_results")
  check(identical(source_data$display_order,1:3) && identical(source_data$colour_class,c("other","TE","agLink")),"SNP_figure_order_and_palette")
  check(settings$value[settings$setting=="bootstrap_replicates"]=="9999" &&
        settings$value[settings$setting=="SNP_density_used"]=="TRUE","SNP_settings_recorded")
  pdf <- file.path(root,"figures/Fig_F_sequential_model_improvement_with_SNP.pdf")
  con <- file(pdf,"rb"); header <- rawToChar(readBin(con,"raw",5)); close(con)
  check(header=="%PDF-" && file.info(pdf)$size>1000,"SNP_figure_is_PDF")
  write_tsv(do.call(rbind,report),file.path(root,"docs/SNP_included_validation_report.tsv"))
  cat("PASS:",length(report),"SNP-included join, model, inference, bootstrap replay and plotting checks.\n")
}


# validate_snp_drop_one
validate_snp_drop_one <- function(root) {
  expected_n <- expected_blocks(root)
  report <- list()
  check <- function(ok,label) {
    report[[length(report)+1L]] <<- data.frame(check=label,status=if(isTRUE(ok)) "PASS" else "FAIL")
    if(!isTRUE(ok)) {
      write_tsv(do.call(rbind,report),file.path(root,"docs/SNP_drop_one_validation_report.tsv"))
      stop("SNP drop-one validation failed: ",label)
    }
  }
  near <- function(a,b,label,tol=1e-9) check(length(a)==length(b) && all(is.finite(a)) &&
    all(is.finite(b)) && max(abs(a-b),0)<tol,label)
  input <- read_tsv(file.path(root,"data/analysis_blocks.tsv.gz"))
  x <- input[input$context=="2Mb", ]; s <- read_tsv(file.path(root,"data/SNP_density_2Mb.tsv.gz"))
  check(nrow(x)==expected_n && nrow(s)==expected_n && !anyDuplicated(x$block_id) && !anyDuplicated(s$block_id),
        "unique_frozen_blocks_match_settings")
  ii <- match(x$block_id,s$block_id)
  check(!anyNA(ii),"SNP_ID_join")
  for(col in c("context","chrom","start","end","context_start","context_end","context_length_bp"))
    check(identical(as.character(x[[col]]),as.character(s[[col]][ii])),paste0("SNP_join_",col))
  z <- function(v) (v-mean(v))/sd(v)
  d <- data.frame(y=log1p(x$length_bp),ag_z=z(log1p(x$aglinks_per_tested_gene)),
    gene_z=z(log1p(x$context_gene_density_per_mb)),recomb_z=z(x$context_recombination_rate),
    te_z=z(x$context_te_content),snp_z=z(log1p(s$context_snp_density_per_mb[ii])),
    acr_z=z(log1p(x$context_acr_density_per_mb)),chrom=factor(x$chrom))
  predictors <- c("ag_z","gene_z","recomb_z","te_z","snp_z","acr_z")
  rows <- read_tsv(file.path(root,"results/model_rows_with_SNP.tsv.gz"))
  check(identical(rows$block_id,x$block_id) && identical(as.character(rows$chrom),as.character(d$chrom)),
        "saved_model_row_IDs_and_chromosomes")
  for(term in c("y",predictors)) near(rows[[term]],d[[term]],paste0("raw_reconstruction_",term))
  X <- model.matrix(~ag_z+gene_z+recomb_z+te_z+snp_z+acr_z+chrom,d)
  full <- lm.fit(X,d$y); full_rss <- sum(full$residuals^2); tss <- sum((d$y-mean(d$y))^2)
  check(full$rank==ncol(X) && nlevels(d$chrom)==21L,"full_rank_and_21_chromosomes")
  tab <- read_tsv(file.path(root,"results/contributions_with_SNP.tsv"))
  p <- read_tsv(file.path(root,"plot_data/Fig_B_drop_one_with_SNP.tsv"))
  nested <- read_tsv(file.path(root,"results/nested_models_with_SNP.tsv"))
  check(nrow(tab)==6L && setequal(tab$removed_term,predictors) && !anyDuplicated(tab$removed_term),
        "six_unique_predictor_removals")
  check(all(tab$n_blocks==expected_n & tab$n_chromosomes==21L & tab$SNP_density_in_full_model &
            tab$chromosome_retained & tab$context=="2Mb" & tab$decomposition=="individual"),
        "identical_analysis_population_and_full_model")
  near(tab$full_R2,rep(1-full_rss/tss,6),"full_R2_from_raw_inputs")
  near(tab$full_R2,rep(nested$R2[nested$stage==3L],6),"same_full_model_as_SNP_sequential_F")
  for(i in seq_len(nrow(tab))) {
    term <- tab$removed_term[i]; label <- paste0(term,": ")
    reduced_X <- X[,colnames(X)!=term,drop=FALSE]
    reduced <- lm.fit(reduced_X,d$y); reduced_rss <- sum(reduced$residuals^2)
    near(c(tab$full_RSS[i],tab$reduced_RSS[i],tab$total_SS[i]),
         c(full_rss,reduced_rss,tss),paste0(label,"residual_and_total_sums_of_squares"),1e-7)
    near(c(tab$full_rank[i],tab$reduced_rank[i]),c(full$rank,reduced$rank),paste0(label,"model_ranks"))
    near(tab$reduced_R2[i],1-reduced_rss/tss,paste0(label,"reduced_R2"))
    near(tab$unique_delta_R2[i],(reduced_rss-full_rss)/tss,paste0(label,"drop_one_delta_R2"))
    near(tab$partial_R2[i],(reduced_rss-full_rss)/reduced_rss,paste0(label,"partial_R2"))
    # A second identity: added variation in y explained by the residualized predictor.
    residual_x <- lm.fit(reduced_X,X[,term])$residuals
    residualized_delta <- as.numeric(crossprod(residual_x,d$y))^2/sum(residual_x^2)/tss
    near(tab$unique_delta_R2[i],residualized_delta,paste0(label,"independent_residualization_identity"))
    chi <- nrow(X)*log(reduced_rss/full_rss)
    near(tab$LRT_chisq[i],chi,paste0(label,"Gaussian_LRT_from_RSS"),1e-8)
    check(tab$LRT_df[i]==1L,paste0(label,"one_removed_parameter"))
    near(log(tab$LRT_P[i]),pchisq(chi,1,lower.tail=FALSE,log.p=TRUE),paste0(label,"LRT_log_P"),1e-7)
    check(setequal(attr(terms(as.formula(tab$full_formula[i])),"term.labels"),c(predictors,"chrom")) &&
          setequal(attr(terms(as.formula(tab$reduced_formula[i])),"term.labels"),c(setdiff(predictors,term),"chrom")),
          paste0(label,"saved_formulas_remove_only_named_predictor"))
  }
  near(tab$unique_delta_R2[tab$removed_term=="ag_z"],nested$step_delta_R2[nested$stage==3L],
       "agLink_drop_one_equals_last_sequential_addition")
  no_snp <- read_tsv(file.path(root,"results/nested_models.tsv"))
  near(tab$reduced_R2[tab$removed_term=="snp_z"],no_snp$R2[no_snp$context=="2Mb" & no_snp$stage==3L],
       "removing_SNP_recovers_existing_SNP_excluded_full_model")
  check(all(tab$unique_delta_R2>=0) && all(diff(tab$unique_delta_R2)<=0),"nonnegative_descending_contributions")
  check(isTRUE(all.equal(p[,names(tab)],tab,tolerance=1e-12,check.attributes=FALSE)),
        "all_plotting_values_match_results")
  check(identical(p$display_order,1:6),"six_bar_display_order")
  near(as.numeric(p$value_label),as.numeric(sprintf("%.4f",tab$unique_delta_R2)),"four_decimal_display_values",1e-12)
  labels <- c(ag_z="agLinks",acr_z="ACR density",gene_z="Gene density",recomb_z="Recombination intensity",
              snp_z="SNP density",te_z="TE content")
  colours <- c(ag_z="#D46A63",acr_z="#E9C46A",gene_z="#999999",recomb_z="#4477AA",snp_z="#999999",te_z="#228833")
  check(identical(p$component,unname(labels[p$removed_term])) &&
        identical(p$colour_hex,unname(colours[p$removed_term])),"predictor_labels_and_reference_palette")
  pdf <- file.path(root,"figures/Fig_B_drop_one_contribution_with_SNP.pdf")
  con <- file(pdf,"rb"); header <- rawToChar(readBin(con,"raw",5)); close(con)
  check(header=="%PDF-" && file.info(pdf)$size>1000,"six_bar_figure_is_PDF")
  write_tsv(do.call(rbind,report),file.path(root,"docs/SNP_drop_one_validation_report.tsv"))
  cat("PASS:",length(report),"SNP drop-one raw-data, independent fit, residualization and figure checks.\n")
}

if(args[1] %in% c("all","plot")) {
  plot_primary(root)
  plot_snp_sequential(root)
  plot_snp_drop_one(root)
}
if(args[1] %in% c("all","validate")) {
  validate_primary(root)
  validate_snp_sequential(root)
  validate_snp_drop_one(root)
}
