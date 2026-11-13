#!/usr/bin/env Rscript
# Multivariate regression from the original frozen block features.
# Consolidates common, fit_models, build_results, refit_models,
# refit_sequential_with_SNP and compute_drop_one_with_SNP.
# Original formulas, transformations, chromosome inference and seeds are retained.
# Shared, base-R-only I/O and figure-data definitions.
read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  read.delim(con, check.names=FALSE, stringsAsFactors=FALSE, na.strings="NA")
}
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  con <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con))
  write.table(x, con, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
input_columns <- c(
  "block_id", "block_name", "chrom", "start", "end", "length_bp",
  "context_start", "context_end", "context_length_bp", "tested_gene_count",
  "aglink_count_gene_block", "aglinks_per_tested_gene", "aglinks_same_block_per_tested_gene",
  "candidate_acr_gene_pair_proxy", "mean_local_acr_per_tested_gene",
  "context_gene_density_per_mb", "context_recombination_rate",
  "context_recomb_coverage_fraction", "context_te_content", "context_acr_density_per_mb"
)
context_order <- c("2Mb", "4Mb")
predictor_terms <- c("ag_z", "gene_z", "recomb_z", "te_z", "acr_z")
contribution_groups <- list(
  individual=list("agLinks"="ag_z", "Gene density"="gene_z",
    "Recombination intensity"="recomb_z", "TE content"="te_z", "ACR density"="acr_z"),
  grouped_C=list("agLinks"="ag_z", "TE content"="te_z",
    "Combined genomic features (ACR + Gene + Recombination)"=c("acr_z", "gene_z", "recomb_z")),
  grouped_E=list("agLinks"="ag_z", "TE content"="te_z", "ACR density"="acr_z",
    "Combined genomic features (Gene + Recombination)"=c("gene_z", "recomb_z"))
)

# This is the sole mapping from statistical results to the four figure files.
make_figure_source <- function(nested, contributions) {
  make_rows <- function(tab, figure, panel, context, labels, colours) {
    n <- nrow(tab)
    out <- data.frame(figure=rep(figure,n), panel=rep(panel,n), context=rep(context,n),
      display_order=seq_len(n), component=labels, colour_class=colours,
      n_blocks=tab$n_blocks, R2=rep(NA_real_,n), unique_delta_R2=rep(NA_real_,n),
      step_delta_R2=rep(NA_real_,n), LRT_P=rep(NA_real_,n),
      wild_bootstrap_P=rep(NA_real_,n), stringsAsFactors=FALSE)
    for (nm in intersect(c("R2","unique_delta_R2","step_delta_R2","LRT_P","wild_bootstrap_P"), names(tab)))
      out[[nm]] <- tab[[nm]]
    out
  }
  models <- nested[nested$context=="2Mb", ]; models <- models[order(models$stage), ]
  stopifnot(nrow(models)==3L)
  A <- make_rows(models[2:3, ], "ABC", "A", "2Mb",
                c("Genomic covariates", "All variables"), c("other","agLink"))
  # A is a two-model comparison: the first bar has no preceding comparison.
  A[1,c("step_delta_R2","LRT_P","wild_bootstrap_P")] <- NA_real_
  rows <- list(A)
  for (spec in list(c("ABC","B","2Mb","individual"), c("ABC","C","2Mb","grouped_C"),
                    c("E","E","2Mb","grouped_E"), c("4Mb","ablation","4Mb","individual"))) {
    tab <- contributions[contributions$context==spec[3] & contributions$decomposition==spec[4], ]
    if (spec[4]=="individual") {
      tab <- tab[order(tab$unique_delta_R2,decreasing=TRUE), ]
      colours <- ifelse(tab$component=="agLinks","agLink",
                   ifelse(tab$component=="TE content","TE",
                     ifelse(tab$component=="Recombination intensity","genomic","other")))
    } else {
      tab <- tab[match(names(contribution_groups[[spec[4]]]),tab$component), ]
      colours <- if(spec[4]=="grouped_C") c("agLink","TE","genomic") else c("agLink","TE","other","genomic")
    }
    rows[[length(rows)+1L]] <- make_rows(tab,spec[1],spec[2],spec[3],tab$component,colours)
  }
  rows[[length(rows)+1L]] <- make_rows(models,"F","F","2Mb",models$model,c("other","TE","agLink"))
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

expected_blocks <- function(root) {
  settings <- read_tsv(file.path(root,"docs/analysis_settings.tsv"))
  n <- as.integer(settings$value[settings$setting=="frozen_blocks_per_context"])
  stopifnot(length(n)==1L,is.finite(n),n>=100L)
  n
}

# Original model inference functions, shared by both analyses.
r2 <- function(m) summary(m)$r.squared
adjr2 <- function(m) summary(m)$adj.r.squared
partial_r2 <- function(m0,m1){
  a=sum(resid(m0)^2); b=sum(resid(m1)^2); (a-b)/a
}
lrt <- function(m0,m1){
  a=logLik(m0); b=logLik(m1)
  chi=2*(as.numeric(b)-as.numeric(a)); df=attr(b,"df")-attr(a,"df")
  c(chisq=chi,df=df,p=pchisq(chi,df,lower.tail=FALSE))
}
safe_inv <- function(A) tryCatch(solve(A),error=function(e) qr.solve(A,diag(nrow(A))))

cluster_vcov <- function(fit,cluster){
  cl=droplevels(factor(cluster)); X=model.matrix(fit); u=resid(fit)
  n=nrow(X); p=ncol(X); G=nlevels(cl); bread=safe_inv(crossprod(X))
  S=matrix(0,ncol(X),G)
  for(g in seq_len(G)){ ii=which(as.integer(cl)==g); S[,g]=crossprod(X[ii,,drop=FALSE],u[ii]) }
  corr=(G/(G-1))*((n-1)/(n-p))
  corr*bread%*%tcrossprod(S)%*%bread
}
cluster_coef <- function(fit,cluster,term){
  V=cluster_vcov(fit,cluster); b=coef(fit); se=sqrt(pmax(0,diag(V)))
  G=nlevels(droplevels(factor(cluster))); tt=b/se
  i=match(term,names(b)); crit=qt(.975,G-1)
  # Strip the coefficient's own name so callers can reliably index beta/se/t/p.
  bi=unname(b[i]); sei=unname(se[i]); ti=unname(tt[i])
  c(beta=bi,se=sei,t=ti,p=unname(2*pt(abs(ti),G-1,lower.tail=FALSE)),
    ci_low=bi-crit*sei,ci_high=bi+crit*sei,G=G)
}
cluster_joint <- function(fit,cluster,terms){
  X=model.matrix(fit); idx=match(terms,colnames(X)); if(anyNA(idx)) stop("Joint terms missing")
  V=cluster_vcov(fit,cluster); b=coef(fit)[idx]; q=length(idx)
  F=as.numeric(crossprod(b,safe_inv(V[idx,idx,drop=FALSE]))%*%b/q)
  G=nlevels(droplevels(factor(cluster)))
  c(F=F,df1=q,df2=G-1,p=pf(F,q,G-1,lower.tail=FALSE))
}
wild_single <- function(full,reduced,cluster,term,B,seed){
  obs=cluster_coef(full,cluster,term)["t"]
  if(B==0) return(list(p=NA_real_,observed=obs,null=numeric()))
  cl=droplevels(factor(cluster)); X=model.matrix(full); qrf=qr(X)
  idx=match(term,colnames(X)); ci=as.integer(cl); G=nlevels(cl)
  fit0=fitted(reduced); u0=resid(reduced)
  # These matrices and row groups are constant over bootstrap replicates.
  n=nrow(X); p=ncol(X); bread=safe_inv(crossprod(X))
  rows=lapply(seq_len(G),function(g) which(ci==g))
  corr=(G/(G-1))*((n-1)/(n-p))
  webb=c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
  out=numeric(B); set.seed(seed)
  for(bi in seq_len(B)){
    w=sample(webb,G,replace=TRUE); ys=fit0+u0*w[ci]
    bs=qr.coef(qrf,ys); us=ys-as.numeric(X%*%bs)
    S=matrix(0,ncol(X),G)
    for(g in seq_len(G)){ ii=rows[[g]]; S[,g]=crossprod(X[ii,,drop=FALSE],us[ii]) }
    V=corr*bread%*%tcrossprod(S)%*%bread
    out[bi]=bs[idx]/sqrt(max(V[idx,idx],0))
  }
  list(p=(1+sum(abs(out)>=abs(obs)))/(B+1),observed=obs,null=out)
}
wild_joint <- function(full,reduced,cluster,terms,B,seed){
  obs=cluster_joint(full,cluster,terms)["F"]
  if(B==0) return(list(p=NA_real_,observed=obs,null=numeric()))
  cl=droplevels(factor(cluster)); X=model.matrix(full); qrf=qr(X)
  idx=match(terms,colnames(X)); ci=as.integer(cl); G=nlevels(cl)
  fit0=fitted(reduced); u0=resid(reduced)
  n=nrow(X); p=ncol(X); bread=safe_inv(crossprod(X))
  rows=lapply(seq_len(G),function(g) which(ci==g))
  corr=(G/(G-1))*((n-1)/(n-p))
  webb=c(-sqrt(1.5),-1,-sqrt(.5),sqrt(.5),1,sqrt(1.5))
  out=numeric(B); set.seed(seed)
  for(bi in seq_len(B)){
    w=sample(webb,G,replace=TRUE); ys=fit0+u0*w[ci]
    bs=qr.coef(qrf,ys); us=ys-as.numeric(X%*%bs)
    S=matrix(0,ncol(X),G)
    for(g in seq_len(G)){ ii=rows[[g]]; S[,g]=crossprod(X[ii,,drop=FALSE],us[ii]) }
    V=corr*bread%*%tcrossprod(S)%*%bread
    bb=bs[idx]; q=length(idx)
    out[bi]=as.numeric(crossprod(bb,safe_inv(V[idx,idx,drop=FALSE]))%*%bb/q)
  }
  list(p=(1+sum(out>=obs))/(B+1),observed=obs,null=out)
}



# fit_context
fit_context <- function(feature_file,outdir,B=9999L,seed=20260910L) {
  RNGkind("Mersenne-Twister","Inversion","Rejection")

  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  tabdir=file.path(outdir,"tables"); plotdir=file.path(outdir,"plot_data")
  figdir=file.path(outdir,"figures"); logdir=file.path(outdir,"logs")
  dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(plotdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(logdir,recursive=TRUE,showWarnings=FALSE)

  write_tsv <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  d <- read.delim(feature_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings=c("NA","NaN",""))

  needed <- c(
   "block_id","chrom","length_bp","tested_gene_count",
   "aglink_count_gene_block","aglinks_per_tested_gene","aglinks_same_block_per_tested_gene",
   "candidate_acr_gene_pair_proxy","mean_local_acr_per_tested_gene",
   "context_gene_density_per_mb","context_recombination_rate","context_recomb_coverage_fraction",
   "context_te_content","context_acr_density_per_mb"
  )
  miss=setdiff(needed,names(d)); if(length(miss)) stop("Missing: ",paste(miss,collapse=", "))
  for(v in setdiff(needed,c("block_id","chrom"))) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))

  # Fixed primary complete-case set.
  keep <- is.finite(d$length_bp) & d$length_bp>0 &
          is.finite(d$tested_gene_count) & d$tested_gene_count>0 &
          is.finite(d$aglinks_per_tested_gene) & d$aglinks_per_tested_gene>=0 &
          is.finite(d$aglinks_same_block_per_tested_gene) & d$aglinks_same_block_per_tested_gene>=0 &
          is.finite(d$context_gene_density_per_mb) & d$context_gene_density_per_mb>=0 &
          is.finite(d$context_recombination_rate) & d$context_recombination_rate>=0 &
          is.finite(d$context_recomb_coverage_fraction) & d$context_recomb_coverage_fraction>=0.5 &
          is.finite(d$context_te_content) & d$context_te_content>=0 & d$context_te_content<=1 &
          is.finite(d$context_acr_density_per_mb) & d$context_acr_density_per_mb>=0 &
          !is.na(d$chrom) & nzchar(d$chrom)

  x=d[keep,,drop=FALSE]
  cnt=table(x$chrom); eligible=names(cnt[cnt>=5])
  x=x[x$chrom %in% eligible,,drop=FALSE]
  x$chrom=droplevels(factor(x$chrom))
  if(nrow(x)<100 || nlevels(x$chrom)<5) stop("Insufficient analysis data")
  if(anyDuplicated(x$block_id)) stop("Duplicate block IDs")

  z <- function(v,nm){
    s=sd(v); if(!is.finite(s)||s<=0) stop("Zero SD: ",nm)
    as.numeric((v-mean(v))/s)
  }

  x$y <- log1p(x$length_bp)
  x$ag_z <- z(log1p(x$aglinks_per_tested_gene),"agLinks/tested gene")
  x$gene_z <- z(log1p(x$context_gene_density_per_mb),"gene density")
  x$recomb_z <- z(x$context_recombination_rate,"recombination")
  x$te_z <- z(x$context_te_content,"TE")
  # Discard any standardized SNP column inherited from a previous frozen table.
  x$snp_z <- NULL
  x$acr_z <- z(log1p(x$context_acr_density_per_mb),"ACR density")

  full_formula <- y ~ ag_z + gene_z + recomb_z + te_z + acr_z + chrom
  reduced_formula <- y ~ gene_z + recomb_z + te_z + acr_z + chrom
  m_full=lm(full_formula,data=x); m_reduced=lm(reduced_formula,data=x)
  if(m_full$rank!=ncol(model.matrix(m_full))) stop("Primary full model rank deficient")

  # Primary inference
  lrA=lrt(m_reduced,m_full)
  cc=cluster_coef(m_full,x$chrom,"ag_z")
  wb=wild_single(m_full,m_reduced,x$chrom,"ag_z",B,seed)

  # A
  A=data.frame(
   model=c("Genomic covariates only","All variables (+ agLinks)"),
   n_blocks=nrow(x),R2=c(r2(m_reduced),r2(m_full)),
   adjusted_R2=c(adjr2(m_reduced),adjr2(m_full)),
   agLink_delta_R2=c(NA,r2(m_full)-r2(m_reduced)),
   agLink_partial_R2=c(NA,partial_r2(m_reduced,m_full)),
   LRT_chisq=c(NA,lrA["chisq"]),LRT_df=c(NA,lrA["df"]),LRT_P=c(NA,lrA["p"]),
   chromosome_cluster_beta=c(NA,cc["beta"]),
   chromosome_cluster_P=c(NA,cc["p"]),
   wild_bootstrap_B=c(NA,B),wild_bootstrap_P=c(NA,wb$p),
   stringsAsFactors=FALSE
  )
  write_tsv(A,file.path(tabdir,"Table2_PanelA_agLink_model_improvement.tsv"))

  # Sequential nested models: baseline -> +TE -> +TE+agLinks. The baseline
  # retains chromosome and every remaining genomic covariate except TE.
  m_baseline=lm(y ~ gene_z + recomb_z + acr_z + chrom,data=x)
  lrTE=lrt(m_baseline,m_reduced)
  te_cc=cluster_coef(m_reduced,x$chrom,"te_z")
  te_wb=wild_single(m_reduced,m_baseline,x$chrom,"te_z",B,seed+4L)
  Ftab=data.frame(
   model=c("baseline","baseline+TE","baseline+TE+agLink"),
   added_predictor=c(NA,"TE content","agLinks"),
   formula=c(paste(deparse(formula(m_baseline)),collapse=" "),
             paste(deparse(formula(m_reduced)),collapse=" "),
             paste(deparse(formula(m_full)),collapse=" ")),
   n_blocks=nrow(x),R2=c(r2(m_baseline),r2(m_reduced),r2(m_full)),
   adjusted_R2=c(adjr2(m_baseline),adjr2(m_reduced),adjr2(m_full)),
   step_delta_R2=c(NA,r2(m_reduced)-r2(m_baseline),r2(m_full)-r2(m_reduced)),
   step_partial_R2=c(NA,partial_r2(m_baseline,m_reduced),partial_r2(m_reduced,m_full)),
   LRT_chisq=c(NA,lrTE["chisq"],lrA["chisq"]),
   LRT_df=c(NA,lrTE["df"],lrA["df"]),LRT_P=c(NA,lrTE["p"],lrA["p"]),
   chromosome_cluster_beta=c(NA,te_cc["beta"],cc["beta"]),
   chromosome_cluster_P=c(NA,te_cc["p"],cc["p"]),
   wild_bootstrap_B=c(NA,B,B),wild_bootstrap_P=c(NA,te_wb$p,wb$p),
   display_order=1:3,colour_class=c("other","TE","agLink")
  )
  write_tsv(Ftab,file.path(tabdir,"Table11_sequential_model_improvement.tsv"))
  write_tsv(Ftab,file.path(plotdir,"Fig_F_sequential_model_improvement.tsv"))
  write_tsv(data.frame(iteration=seq_len(B),null_cluster_t=te_wb$null),
            file.path(tabdir,"Table12_TE_addition_wild_null.tsv"))

  # Full coefficients (standardized predictors, unstandardized log-length response)
  ct=coef(summary(m_full))
  terms=c(ag_z="agLinks",gene_z="Gene density",recomb_z="Recombination intensity",
          te_z="TE content",acr_z="ACR density")
  coefout=do.call(rbind,lapply(names(terms),function(tn){
    cr=cluster_coef(m_full,x$chrom,tn)
    data.frame(predictor=unname(terms[tn]),term=tn,beta=ct[tn,"Estimate"],
               OLS_SE=ct[tn,"Std. Error"],OLS_P=ct[tn,"Pr(>|t|)"],
               chromosome_cluster_SE=cr["se"],chromosome_cluster_P=cr["p"],
               cluster_CI_low=cr["ci_low"],cluster_CI_high=cr["ci_high"])
  }))
  write_tsv(coefout,file.path(tabdir,"Table1_full_model_coefficients.tsv"))

  # B drop-one
  drop_groups=list(
   "agLinks"="ag_z",
   "Gene density"="gene_z",
   "Recombination intensity"="recomb_z",
   "ACR density"="acr_z",
   "TE content"="te_z"
  )
  all_terms=c("ag_z","gene_z","recomb_z","te_z","acr_z")
  Btab=do.call(rbind,lapply(names(drop_groups),function(nm){
    kept=setdiff(all_terms,drop_groups[[nm]])
    mr=lm(reformulate(c(kept,"chrom"),response="y"),data=x)
    ll=lrt(mr,m_full)
    data.frame(predictor=nm,n_blocks=nrow(x),full_R2=r2(m_full),reduced_R2=r2(mr),
               unique_delta_R2=r2(m_full)-r2(mr),partial_R2=partial_r2(mr,m_full),
               LRT_chisq=ll["chisq"],LRT_df=ll["df"],LRT_P=ll["p"])
  }))
  Btab=Btab[order(Btab$unique_delta_R2,decreasing=TRUE),]
  write_tsv(Btab,file.path(tabdir,"Table3_PanelB_drop_one_contributions.tsv"))

  # C grouped drop-one
  Cgroups=list(
   "agLinks"="ag_z",
   "TE content"="te_z",
   "Combined genomic features (ACR + Gene + Recombination)"=
     c("acr_z","gene_z","recomb_z")
  )
  Ctab=do.call(rbind,lapply(names(Cgroups),function(nm){
    kept=setdiff(all_terms,Cgroups[[nm]])
    mr=lm(reformulate(c(kept,"chrom"),response="y"),data=x)
    ll=lrt(mr,m_full)
    data.frame(component=nm,n_blocks=nrow(x),full_R2=r2(m_full),reduced_R2=r2(mr),
               unique_delta_R2=r2(m_full)-r2(mr),partial_R2=partial_r2(mr,m_full),
               LRT_chisq=ll["chisq"],LRT_df=ll["df"],LRT_P=ll["p"])
  }))
  write_tsv(Ctab,file.path(tabdir,"Table4_PanelC_grouped_contributions.tsv"))

  # E grouped drop-one: separate ACR density from the other genomic features.
  Egroups=list(
   "agLinks"="ag_z",
   "TE content"="te_z",
   "ACR density"="acr_z",
   "Combined genomic features (Gene + Recombination)"=
     c("gene_z","recomb_z")
  )
  Etab=do.call(rbind,lapply(names(Egroups),function(nm){
    kept=setdiff(all_terms,Egroups[[nm]])
    mr=lm(reformulate(c(kept,"chrom"),response="y"),data=x)
    ll=lrt(mr,m_full)
    data.frame(component=nm,n_blocks=nrow(x),full_R2=r2(m_full),reduced_R2=r2(mr),
               unique_delta_R2=r2(m_full)-r2(mr),partial_R2=partial_r2(mr,m_full),
               LRT_chisq=ll["chisq"],LRT_df=ll["df"],LRT_P=ll["p"])
  }))
  write_tsv(Etab,file.path(tabdir,"Table10_PanelE_grouped_contributions.tsv"))

  # Sensitivity 1: two-part, same covariates, NO extra opportunity splines.
  x$has_aglink=as.integer(x$aglinks_per_tested_gene>0)
  pos=x$has_aglink==1
  x$positive_log_degree_centered=0
  center=median(log(x$aglinks_per_tested_gene[pos]))
  x$positive_log_degree_centered[pos]=log(x$aglinks_per_tested_gene[pos])-center
  two_terms=c("has_aglink","positive_log_degree_centered")
  m2r=lm(y ~ gene_z+recomb_z+te_z+acr_z+chrom,data=x)
  m2f=lm(y ~ has_aglink+positive_log_degree_centered+gene_z+recomb_z+te_z+acr_z+chrom,data=x)
  j2=cluster_joint(m2f,x$chrom,two_terms)
  w2=wild_joint(m2f,m2r,x$chrom,two_terms,B,seed+1L)

  # Sensitivity 2: same-block degree, same primary covariates.
  x$ag_same_z=z(log1p(x$aglinks_same_block_per_tested_gene),"same-block agLink degree")
  msr=m_reduced
  msf=lm(y ~ ag_same_z+gene_z+recomb_z+te_z+acr_z+chrom,data=x)
  cs=cluster_coef(msf,x$chrom,"ag_same_z")
  ws=wild_single(msf,msr,x$chrom,"ag_same_z",B,seed+2L)

  # Sensitivity 3: stringent opportunity-aware two-part model.
  # Uses only rows from the primary analysis set that have positive candidate-pair opportunity.
  xs=x[is.finite(x$candidate_acr_gene_pair_proxy) & x$candidate_acr_gene_pair_proxy>0 &
       is.finite(x$tested_gene_count) & x$tested_gene_count>0,,drop=FALSE]
  xs$log_tested=log(xs$tested_gene_count)
  xs$log_pairs=log(xs$candidate_acr_gene_pair_proxy)
  xs$has_aglink=as.integer(xs$aglinks_per_tested_gene>0)
  pp=xs$has_aglink==1
  xs$positive_log_degree_centered=0
  ccen=median(log(xs$aglinks_per_tested_gene[pp]))
  xs$positive_log_degree_centered[pp]=log(xs$aglinks_per_tested_gene[pp])-ccen
  m3r=lm(y ~ splines::ns(log_tested,df=4)+splines::ns(log_pairs,df=3)+
            gene_z+recomb_z+te_z+acr_z+chrom,data=xs)
  m3f=lm(y ~ has_aglink+positive_log_degree_centered+
            splines::ns(log_tested,df=4)+splines::ns(log_pairs,df=3)+
            gene_z+recomb_z+te_z+acr_z+chrom,data=xs)
  j3=cluster_joint(m3f,xs$chrom,two_terms)
  w3=wild_joint(m3f,m3r,xs$chrom,two_terms,B,seed+3L)

  Sens=data.frame(
   analysis=c(
     "PRIMARY: log1p(agLinks/tested gene), genomic covariates + chromosome",
     "Two-part agLink, same genomic covariates + chromosome",
     "Same-block agLinks/tested gene, same genomic covariates + chromosome",
     "STRINGENT: two-part + tested-gene/candidate-pair opportunity splines"
   ),
   n_blocks=c(nrow(x),nrow(x),nrow(x),nrow(xs)),
   delta_R2=c(
     r2(m_full)-r2(m_reduced),
     r2(m2f)-r2(m2r),
     r2(msf)-r2(msr),
     r2(m3f)-r2(m3r)
   ),
   partial_R2=c(
     partial_r2(m_reduced,m_full),
     partial_r2(m2r,m2f),
     partial_r2(msr,msf),
     partial_r2(m3r,m3f)
   ),
   chromosome_cluster_P=c(cc["p"],j2["p"],cs["p"],j3["p"]),
   wild_cluster_P=c(wb$p,w2$p,ws$p,w3$p),
   stringsAsFactors=FALSE
  )
  write_tsv(Sens,file.path(tabdir,"Table5_agLink_sensitivity_hierarchy.tsv"))

  # Correlation + VIF
  corrdat=data.frame(
   block_length=x$length_bp,aglinks_per_tested_gene=x$aglinks_per_tested_gene,
   gene_density=x$context_gene_density_per_mb,recombination=x$context_recombination_rate,
   TE=x$context_te_content,
   ACR_density=x$context_acr_density_per_mb,tested_gene_count=x$tested_gene_count,
   candidate_pair_proxy=x$candidate_acr_gene_pair_proxy
  )
  cm=cor(corrdat,method="spearman",use="pairwise.complete.obs")
  write_tsv(data.frame(variable=rownames(cm),cm,check.names=FALSE),
            file.path(tabdir,"Table6_Spearman_correlation.tsv"))

  vifs=data.frame(term=all_terms,VIF=NA_real_)
  for(i in seq_along(all_terms)){
    response=all_terms[i]; others=c(setdiff(all_terms,response),"chrom")
    mm=lm(reformulate(others,response=response),data=x)
    vifs$VIF[i]=1/(1-summary(mm)$r.squared)
  }
  write_tsv(vifs,file.path(tabdir,"Table7_VIF.tsv"))
  write_tsv(x,file.path(tabdir,"Table0_complete_case_blocks.tsv"))
  write_tsv(data.frame(iteration=seq_len(B),null_cluster_t=wb$null),
            file.path(tabdir,"Table8_primary_agLink_wild_null.tsv"))

  # Plot data
  Ap=A; Ap$display_label=c("Genomic covariates","All variables"); Ap$display_order=1:2
  write_tsv(Ap,file.path(plotdir,"Fig_ABC_panel_A_model_improvement.tsv"))
  Bp=Btab; Bp$display_order=seq_len(nrow(Bp))
  Bp$colour_class=ifelse(Bp$predictor=="agLinks","agLink",
                   ifelse(Bp$predictor=="Recombination intensity","genomic",
                   ifelse(Bp$predictor=="TE content","TE","other")))
  write_tsv(Bp,file.path(plotdir,"Fig_ABC_panel_B_drop_one.tsv"))
  ord=c("agLinks","TE content","Combined genomic features (ACR + Gene + Recombination)")
  Cp=Ctab[match(ord,Ctab$component),]; Cp$display_order=1:3
  Cp$colour_class=c("agLink","TE","genomic")
  write_tsv(Cp,file.path(plotdir,"Fig_ABC_panel_C_grouped.tsv"))
  ordE=c("agLinks","TE content","ACR density",
         "Combined genomic features (Gene + Recombination)")
  Ep=Etab[match(ordE,Etab$component),]; Ep$display_order=1:4
  Ep$colour_class=c("agLink","TE","other","genomic")
  write_tsv(Ep,file.path(plotdir,"Fig_E_panel_grouped.tsv"))

  # QC + summary
  QC=data.frame(
   metric=c("n_blocks","n_chromosomes","primary_reduced_R2","primary_full_R2",
            "primary_agLink_delta_R2","primary_agLink_partial_R2",
            "primary_LRT_P","primary_chr_cluster_P","primary_wild_cluster_P",
            "two_part_degree_center"),
   value=c(nrow(x),nlevels(x$chrom),r2(m_reduced),r2(m_full),
           r2(m_full)-r2(m_reduced),partial_r2(m_reduced,m_full),
           lrA["p"],cc["p"],wb$p,exp(center))
  )
  write_tsv(QC,file.path(tabdir,"Table9_analysis_QC.tsv"))

  fmtp=function(p) if(!is.finite(p)) "NA" else if(p<1e-4) format(p,scientific=TRUE,digits=3) else sprintf("%.4f",p)
  lines=c(
   "FINAL REVIEWER-FACING MULTIVARIATE REGRESSION",
   paste0("Blocks: ",nrow(x),"; chromosomes: ",nlevels(x$chrom)),
   "Primary response: log1p(haplotype-block length).",
   "Primary agLink predictor: z[log1p(unique significant agLinks per tested gene)].",
   "Genomic covariates: fixed-context gene density, recombination consensus rank, TE fraction, ACR density. SNP density excluded.",
   "Chromosome is a categorical fixed effect.",
   "Sequential baseline: gene density + recombination consensus rank + ACR density + chromosome.",
   "",
   paste0("Reduced R2=",sprintf("%.4f",r2(m_reduced)),
          "; full R2=",sprintf("%.4f",r2(m_full)),
          "; agLink unique delta R2=",sprintf("%.4f",r2(m_full)-r2(m_reduced)),
          "; partial R2=",sprintf("%.4f",partial_r2(m_reduced,m_full))),
   paste0("Ordinary nested-model LRT P=",fmtp(lrA["p"])),
   paste0("Chromosome-clustered P=",fmtp(cc["p"])),
   paste0("Webb wild chromosome-cluster P=",fmtp(wb$p)," (B=",B,")"),
   paste0("Sequential R2: baseline=",sprintf("%.4f",r2(m_baseline)),
          "; baseline+TE=",sprintf("%.4f",r2(m_reduced)),
          "; baseline+TE+agLink=",sprintf("%.4f",r2(m_full))),
   paste0("Adding TE to baseline: delta R2=",sprintf("%.4f",r2(m_reduced)-r2(m_baseline)),
          "; LRT P=",fmtp(lrTE["p"]),"; wild P=",fmtp(te_wb$p)),
   "",
   "Interpretation: adjusted association beyond the included genomic covariates; not proof of causality or selection.",
   "The stringent opportunity-aware two-part model is retained as sensitivity rather than used as the primary reviewer-facing estimand."
  )
  writeLines(lines,file.path(outdir,"reviewer_ready_summary.txt"))
  capture.output(sessionInfo(),file=file.path(logdir,"R_sessionInfo.txt"))
  cat("DONE\n")
}

# Combine two temporary model-fit directories into the canonical result tables.
# Loaded by refit_models.R; it can also import a previously verified fit.
assemble_results <- function(fit_dirs, root, seed) {
  stopifnot(identical(names(fit_dirs), context_order))
  get_table <- function(context, filename) read_tsv(file.path(fit_dirs[[context]],"tables",filename))
  add_context <- function(tab, context) data.frame(context=context,tab,check.names=FALSE,stringsAsFactors=FALSE)
  stack <- function(filename) do.call(rbind,lapply(context_order,function(ctx) add_context(get_table(ctx,filename),ctx)))

  nested <- stack("Table11_sequential_model_improvement.tsv")
  nested$stage <- nested$display_order; nested$display_order <- NULL; nested$colour_class <- NULL
  nested <- nested[,c("context","stage",setdiff(names(nested),c("context","stage")))]
  coefficients <- stack("Table1_full_model_coefficients.tsv")
  sensitivity <- stack("Table5_agLink_sensitivity_hierarchy.tsv")
  contributions <- list(); diagnostics <- list(); nulls <- list()
  for (ctx in context_order) {
    for (spec in list(c("individual","Table3_PanelB_drop_one_contributions.tsv"),
                      c("grouped_C","Table4_PanelC_grouped_contributions.tsv"),
                      c("grouped_E","Table10_PanelE_grouped_contributions.tsv"))) {
      tab <- get_table(ctx,spec[2]); names(tab)[names(tab)=="predictor"] <- "component"
      contributions[[length(contributions)+1L]] <- data.frame(context=ctx,decomposition=spec[1],tab,
                                                              check.names=FALSE,stringsAsFactors=FALSE)
    }
    qc <- get_table(ctx,"Table9_analysis_QC.tsv")
    diagnostics[[length(diagnostics)+1L]] <- data.frame(context=ctx,diagnostic="analysis_QC",
      variable_1=qc$metric,variable_2=NA_character_,value=qc$value,stringsAsFactors=FALSE)
    vif <- get_table(ctx,"Table7_VIF.tsv")
    diagnostics[[length(diagnostics)+1L]] <- data.frame(context=ctx,diagnostic="VIF",
      variable_1=vif$term,variable_2=NA_character_,value=vif$VIF,stringsAsFactors=FALSE)
    corr <- get_table(ctx,"Table6_Spearman_correlation.tsv")
    # Retain each correlation once; omit redundant diagonal and mirror entries.
    mat <- as.matrix(corr[,-1]); ij <- which(upper.tri(mat),arr.ind=TRUE)
    diagnostics[[length(diagnostics)+1L]] <- data.frame(context=ctx,diagnostic="Spearman",
      variable_1=corr$variable[ij[,1]],variable_2=names(corr)[-1][ij[,2]],value=mat[ij],stringsAsFactors=FALSE)
    for (spec in list(c("add_agLink","Table8_primary_agLink_wild_null.tsv"),
                      c("add_TE","Table12_TE_addition_wild_null.tsv"))) {
      tab <- get_table(ctx,spec[2])
      nulls[[length(nulls)+1L]] <- data.frame(context=ctx,comparison=spec[1],tab,stringsAsFactors=FALSE)
    }
  }
  contributions <- do.call(rbind,contributions)
  write_tsv(coefficients,file.path(root,"results/model_coefficients.tsv"))
  write_tsv(nested,file.path(root,"results/nested_models.tsv"))
  write_tsv(contributions,file.path(root,"results/contributions.tsv"))
  write_tsv(sensitivity,file.path(root,"results/sensitivity.tsv"))
  write_tsv(do.call(rbind,diagnostics),file.path(root,"results/diagnostics.tsv"))
  write_tsv(do.call(rbind,nulls),file.path(root,"results/bootstrap_null.tsv.gz"))
  write_tsv(make_figure_source(nested,contributions),file.path(root,"plot_data/figure_source_data.tsv"))
  settings <- data.frame(setting=c("model_version","contexts","bootstrap_replicates","random_seed",
    "RNG_uniform","RNG_normal","RNG_sample","SNP_density_used","frozen_blocks_per_context"),
    value=c("no_SNP_20261019",paste(context_order,collapse=","),
      unique(nested$wild_bootstrap_B[is.finite(nested$wild_bootstrap_B)]),seed,
      "Mersenne-Twister","Inversion","Rejection","FALSE",unique(nested$n_blocks)),stringsAsFactors=FALSE)
  write_tsv(settings,file.path(root,"docs/analysis_settings.tsv"))
  sessions <- unlist(lapply(context_order,function(ctx) c(paste0("CONTEXT: ",ctx),
    readLines(file.path(fit_dirs[[ctx]],"logs/R_sessionInfo.txt")),"")))
  writeLines(sessions,file.path(root,"docs/R_sessionInfo.txt"))
  summary_lines <- c("FINAL MULTIVARIATE REGRESSION (SNP DENSITY EXCLUDED)",
    sprintf("Frozen block-level features; %d blocks per context, 21 chromosomes.",unique(nested$n_blocks)),
    "R2 is in-sample; all effects are adjusted associations, not proof of causality.", "")
  for(ctx in context_order) {
    nm <- nested[nested$context==ctx, ]; nm <- nm[order(nm$stage), ]
    st <- sensitivity[sensitivity$context==ctx, ]
    summary_lines <- c(summary_lines,paste("CONTEXT:",ctx),
      paste0("R2: baseline=",sprintf("%.7f",nm$R2[1]),"; baseline+TE=",sprintf("%.7f",nm$R2[2]),
             "; baseline+TE+agLink=",sprintf("%.7f",nm$R2[3])),
      paste0("TE addition: delta R2=",sprintf("%.7f",nm$step_delta_R2[2]),
             "; LRT P=",format(nm$LRT_P[2],digits=5),"; wild P=",sprintf("%.4f",nm$wild_bootstrap_P[2])),
      paste0("agLink addition: delta R2=",sprintf("%.7f",nm$step_delta_R2[3]),
             "; LRT P=",format(nm$LRT_P[3],digits=5),"; wild P=",sprintf("%.4f",nm$wild_bootstrap_P[3])),
      paste0("Stringent opportunity-adjusted sensitivity: delta R2=",sprintf("%.7f",st$delta_R2[4]),
             "; wild P=",sprintf("%.4f",st$wild_cluster_P[4])),"")
  }
  writeLines(summary_lines,file.path(root,"docs/RESULTS_SUMMARY.txt"))
  invisible(TRUE)
}


# refit_primary
refit_primary <- function(root) {
  settings <- read_tsv(file.path(root,"docs/analysis_settings.tsv"))
  get_setting <- function(key) settings$value[match(key,settings$setting)]
  B <- as.integer(get_setting("bootstrap_replicates"))
  seed <- as.integer(get_setting("random_seed"))
  stopifnot(B==9999L,is.finite(seed))
  input <- read_tsv(file.path(root,"data/analysis_blocks.tsv.gz"))
  stopifnot(setequal(unique(input$context),context_order))
  tmp <- tempfile("multireg_submission_refit_"); dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  fit_dirs <- setNames(file.path(tmp,context_order),context_order)
  for (ctx in context_order) {
    cat("Refitting",ctx,"with",B,"wild-bootstrap replicates...\n")
    feature_file <- file.path(tmp,paste0(ctx,"_input.tsv"))
    write_tsv(input[input$context==ctx,input_columns],feature_file)
    fit_context(feature_file,fit_dirs[[ctx]],B,seed)
  }
  # Intermediate legacy-shaped tables are temporary; only consolidated tables are published.
  assemble_results(fit_dirs,root,seed)
  unlink(tmp,recursive=TRUE)
  cat("DONE: both contexts refitted; statistical and figure-source tables regenerated.\n")
}


# refit_with_snp
refit_with_snp <- function(root,B=9999L) {
  out <- root
  expected_n <- expected_blocks(root)
  seed <- 20260910L
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  input <- read_tsv(file.path(root,"data/analysis_blocks.tsv.gz"))
  x <- input[input$context=="2Mb", ]; snp <- read_tsv(file.path(root,"data/SNP_density_2Mb.tsv.gz"))
  stopifnot(nrow(x)==expected_n,nrow(snp)==expected_n,!anyDuplicated(snp$block_id))
  ii <- match(x$block_id,snp$block_id); stopifnot(!anyNA(ii))
  for(col in c("context","block_id","chrom","start","end","context_start","context_end","context_length_bp"))
    stopifnot(identical(as.character(x[[col]]),as.character(snp[[col]][ii])))
  x$context_snp_density_per_mb <- snp$context_snp_density_per_mb[ii]
  x$context_snp_count <- snp$context_snp_count[ii]
  stopifnot(all(is.finite(x$context_snp_density_per_mb)),all(x$context_snp_density_per_mb>=0))
  z <- function(v) { stopifnot(is.finite(sd(v)),sd(v)>0); as.numeric(scale(v)) }
  x$chrom <- droplevels(factor(x$chrom)); stopifnot(nlevels(x$chrom)==21L)
  x$y <- log1p(x$length_bp)
  x$ag_z <- z(log1p(x$aglinks_per_tested_gene))
  x$gene_z <- z(log1p(x$context_gene_density_per_mb))
  x$recomb_z <- z(x$context_recombination_rate)
  x$acr_z <- z(log1p(x$context_acr_density_per_mb))
  x$snp_z <- z(log1p(x$context_snp_density_per_mb))
  x$te_z <- z(x$context_te_content)
  models <- list(
    lm(y~gene_z+recomb_z+snp_z+acr_z+chrom,data=x),
    lm(y~gene_z+recomb_z+te_z+snp_z+acr_z+chrom,data=x),
    lm(y~ag_z+gene_z+recomb_z+te_z+snp_z+acr_z+chrom,data=x))
  stopifnot(all(vapply(models,function(m) m$rank==ncol(model.matrix(m)),logical(1))))
  labels <- c("baseline","baseline+TE","baseline+TE+agLink")
  tab <- data.frame(context="2Mb",stage=1:3,model=labels,
    added_predictor=c(NA,"TE content","agLinks"),
    formula=vapply(models,function(m) paste(deparse(formula(m)),collapse=" "),character(1)),
    n_blocks=nrow(x),n_chromosomes=nlevels(x$chrom),SNP_density_used=TRUE,
    R2=vapply(models,r2,numeric(1)),adjusted_R2=vapply(models,adjr2,numeric(1)),
    step_delta_R2=NA_real_,step_partial_R2=NA_real_,LRT_chisq=NA_real_,LRT_df=NA_integer_,LRT_P=NA_real_,
    chromosome_cluster_beta=NA_real_,chromosome_cluster_P=NA_real_,
    observed_cluster_t=NA_real_,wild_bootstrap_B=NA_integer_,wild_bootstrap_seed=NA_integer_,
    wild_bootstrap_P=NA_real_,stringsAsFactors=FALSE)
  nulls <- list()
  for(k in 2:3) {
    term <- if(k==2) "te_z" else "ag_z"
    current_seed <- if(k==2) seed+4L else seed
    ll <- lrt(models[[k-1]],models[[k]])
    cc <- cluster_coef(models[[k]],x$chrom,term)
    cat("With SNP: testing",labels[k-1],"->",labels[k],"with",B,"Webb replicates...\n")
    wb <- wild_single(models[[k]],models[[k-1]],x$chrom,term,B,current_seed)
    tab$step_delta_R2[k] <- tab$R2[k]-tab$R2[k-1]
    tab$step_partial_R2[k] <- partial_r2(models[[k-1]],models[[k]])
    tab[k,c("LRT_chisq","LRT_df","LRT_P")] <- unname(ll)
    tab$chromosome_cluster_beta[k] <- cc["beta"]
    tab$chromosome_cluster_P[k] <- cc["p"]
    tab$observed_cluster_t[k] <- wb$observed
    tab$wild_bootstrap_B[k] <- B; tab$wild_bootstrap_seed[k] <- current_seed
    tab$wild_bootstrap_P[k] <- wb$p
    nulls[[length(nulls)+1L]] <- data.frame(context="2Mb",SNP_density_used=TRUE,
      comparison=if(k==2) "add_TE" else "add_agLink",iteration=seq_len(B),null_cluster_t=wb$null)
  }
  coeff <- do.call(rbind,lapply(seq_along(models),function(k) {
    m <- models[[k]]; ct <- coef(summary(m)); V <- cluster_vcov(m,x$chrom)
    se <- sqrt(pmax(0,diag(V))); crit <- qt(.975,nlevels(x$chrom)-1)
    answer <- data.frame(stage=k,model=labels[k],term=rownames(ct),beta=unname(ct[,1]),
      OLS_SE=unname(ct[,2]),OLS_P=unname(ct[,4]),chromosome_cluster_SE=unname(se),
      chromosome_cluster_P=unname(2*pt(abs(ct[,1]/se),nlevels(x$chrom)-1,lower.tail=FALSE)),
      cluster_CI_low=unname(ct[,1]-crit*se),cluster_CI_high=unname(ct[,1]+crit*se),row.names=NULL)
    # Match the primary compact tables: chromosome/intercept are fitted but
    # omitted from the predictor coefficient table (their formulas remain saved).
    answer[answer$term %in% c("ag_z","gene_z","recomb_z","te_z","snp_z","acr_z"), ]
  }))
  write_tsv(tab,file.path(out,"results/nested_models_with_SNP.tsv"))
  write_tsv(coeff,file.path(out,"results/model_coefficients_with_SNP.tsv"))
  write_tsv(do.call(rbind,nulls),file.path(out,"results/bootstrap_null_with_SNP.tsv.gz"))
  row_cols <- c("block_id","chrom","start","end","length_bp","context_start","context_end","context_length_bp",
    "aglinks_per_tested_gene","context_gene_density_per_mb","context_recombination_rate",
    "context_te_content","context_acr_density_per_mb","context_snp_count","context_snp_density_per_mb",
    "y","ag_z","gene_z","recomb_z","te_z","acr_z","snp_z")
  write_tsv(x[,row_cols],file.path(out,"results/model_rows_with_SNP.tsv.gz"))
  plot <- tab; plot$display_order <- plot$stage; plot$colour_class <- c("other","TE","agLink")
  write_tsv(plot,file.path(out,"plot_data/Fig_F_sequential_with_SNP.tsv"))
  settings <- data.frame(setting=c("model_version","context","bootstrap_replicates","random_seed",
    "TE_addition_seed","agLink_addition_seed","SNP_density_used","SNP_transformation",
    "n_blocks","n_chromosomes","RNG_uniform","RNG_normal","RNG_sample"),
    value=c("additional_F_with_SNP_20261027","2Mb",B,seed,seed+4L,seed,"TRUE",
      "z(log1p(context_snp_density_per_mb))",nrow(x),nlevels(x$chrom),
      "Mersenne-Twister","Inversion","Rejection"))
  write_tsv(settings,file.path(out,"docs/SNP_included_settings.tsv"))
  capture.output(sessionInfo(),file=file.path(out,"docs/SNP_included_R_sessionInfo.txt"))
  lines <- c("ADDITIONAL F FIGURE: SNP DENSITY INCLUDED IN EVERY MODEL",
    sprintf("Original 2-Mb inputs: %d blocks, 21 chromosomes.",nrow(x)),
    "Baseline: Gene density + Recombination + ACR density + SNP density + chromosome.",
    paste("R2:",paste(sprintf("%s=%.7f",tab$model,tab$R2),collapse="; ")),
    sprintf("TE addition: delta R2=%.7f; LRT P=%.8g; wild P=%.4f",tab$step_delta_R2[2],tab$LRT_P[2],tab$wild_bootstrap_P[2]),
    sprintf("agLink addition: delta R2=%.7f; LRT P=%.8g; wild P=%.4f",tab$step_delta_R2[3],tab$LRT_P[3],tab$wild_bootstrap_P[3]),
    "This is an additional SNP-adjusted analysis; the SNP-excluded primary figures/results are retained.",
    "R2 is ordinary in-sample fit. Wild tests use chromosome-level Webb weights, not CRU sign flips.")
  writeLines(lines,file.path(out,"docs/SNP_included_RESULTS_SUMMARY.txt"))
  cat(paste(lines,collapse="\n"),"\n")
}


# compute_snp_drop_one
compute_snp_drop_one <- function(root) {
  out <- root
  expected_n <- expected_blocks(root)
  x <- read_tsv(file.path(root,"results/model_rows_with_SNP.tsv.gz"))
  nested <- read_tsv(file.path(root,"results/nested_models_with_SNP.tsv"))
  x$chrom <- factor(x$chrom)
  terms <- c("ag_z","gene_z","recomb_z","te_z","snp_z","acr_z")
  stopifnot(nrow(x)==expected_n,!anyDuplicated(x$block_id),nlevels(x$chrom)==21L,
            all(is.finite(as.matrix(x[,c("y",terms)]))))
  full <- lm(y~ag_z+gene_z+recomb_z+te_z+snp_z+acr_z+chrom,data=x,na.action=na.fail)
  stopifnot(full$rank==ncol(model.matrix(full)),
            abs(summary(full)$r.squared-nested$R2[nested$stage==3L])<1e-12)
  labels <- c(ag_z="agLinks",acr_z="ACR density",gene_z="Gene density",
              recomb_z="Recombination intensity",snp_z="SNP density",te_z="TE content")
  colours <- c(ag_z="#D46A63",acr_z="#E9C46A",gene_z="#999999",
               recomb_z="#4477AA",snp_z="#999999",te_z="#228833")
  full_rss <- sum(resid(full)^2); total_ss <- sum((x$y-mean(x$y))^2)
  tab <- do.call(rbind,lapply(names(labels),function(term) {
    reduced <- lm(reformulate(c(setdiff(terms,term),"chrom"),response="y"),data=x,na.action=na.fail)
    stopifnot(nobs(reduced)==nobs(full),reduced$rank==full$rank-1L,
              identical(rownames(model.frame(full)),rownames(model.frame(reduced))))
    reduced_rss <- sum(resid(reduced)^2)
    chi <- 2*(as.numeric(logLik(full))-as.numeric(logLik(reduced)))
    df <- attr(logLik(full),"df")-attr(logLik(reduced),"df")
    data.frame(context="2Mb",decomposition="individual",component=unname(labels[term]),
      removed_term=term,n_blocks=nrow(x),n_chromosomes=nlevels(x$chrom),
      SNP_density_in_full_model=TRUE,chromosome_retained=TRUE,
      full_formula=paste(deparse(formula(full)),collapse=" "),
      reduced_formula=paste(deparse(formula(reduced)),collapse=" "),
      full_rank=full$rank,reduced_rank=reduced$rank,
      full_RSS=full_rss,reduced_RSS=reduced_rss,total_SS=total_ss,
      full_R2=summary(full)$r.squared,reduced_R2=summary(reduced)$r.squared,
      unique_delta_R2=summary(full)$r.squared-summary(reduced)$r.squared,
      partial_R2=1-full_rss/reduced_rss,LRT_chisq=chi,LRT_df=df,
      LRT_P=pchisq(chi,df,lower.tail=FALSE),stringsAsFactors=FALSE)
  }))
  tab <- tab[order(tab$unique_delta_R2,decreasing=TRUE), ]; rownames(tab) <- NULL
  stopifnot(nrow(tab)==6L,all(tab$unique_delta_R2>=0))
  write_tsv(tab,file.path(out,"results/contributions_with_SNP.tsv"))
  plot <- tab
  plot$display_order <- seq_len(nrow(plot))
  plot$colour_hex <- unname(colours[plot$removed_term])
  plot$value_label <- sprintf("%.4f",plot$unique_delta_R2)
  write_tsv(plot,file.path(out,"plot_data/Fig_B_drop_one_with_SNP.tsv"))
  lines <- c("SNP-INCLUDED 2-Mb INDIVIDUAL DROP-ONE CONTRIBUTIONS",
    sprintf("Same %d blocks and 21 chromosomes as the SNP-included sequential F figure.",nrow(x)),
    paste("Full model:",tab$full_formula[1]),sprintf("Full R2: %.12f",tab$full_R2[1]),
    "Each reduced model removes exactly one predictor; chromosome is always retained.",
    "Unique delta R2 = full-model R2 - reduced-model R2 (not sequential addition).",
    sprintf("%s: unique delta R2 = %.12f; plotted label = %.4f",tab$component,tab$unique_delta_R2,tab$unique_delta_R2),
    "LRT_P in the result table is the ordinary Gaussian-model LRT, not a wild-bootstrap P.",
    "The figure displays effect sizes only. Contributions are not additive or causal.")
  writeLines(lines,file.path(out,"docs/SNP_drop_one_RESULTS_SUMMARY.txt"))
  capture.output(sessionInfo(),file=file.path(out,"docs/SNP_drop_one_R_sessionInfo.txt"))
  cat(paste(lines,collapse="\n"),"\n")
}

run_models <- function(root) {
  refit_primary(root)
  refit_with_snp(root)
  compute_snp_drop_one(root)
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=1L) stop("Usage: Rscript 01_fit_models.R OUTPUT_DIR")
  run_models(normalizePath(args[1],mustWork=TRUE))
}
