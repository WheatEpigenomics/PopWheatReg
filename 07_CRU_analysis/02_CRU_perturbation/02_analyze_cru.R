#!/usr/bin/env Rscript
# CRU perturbation: statistics, figures and validation (base R only).
# Consolidated from compute/plot/validate_sentinel_contrasts.R.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L || !args[1] %in% c("all","statistics","plot","validate")) {
  stop("Usage: Rscript 02_analyze_cru.R {all|statistics|plot|validate} OUTPUT_DIR")
}
root <- normalizePath(args[2],mustWork=TRUE)
rd <- function(name) {
  path <- file.path(root,name)
  con <- if(grepl("[.]gz$",path)) gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con))
  read.delim(con,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA","nan"))
}
wr <- function(x,name) {
  path <- file.path(root,name)
  con <- if(grepl("[.]gz$",path)) gzfile(path,"wt") else file(path,"wt")
  on.exit(close(con)); write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}


# ========================================================================
# compute_statistics
# ========================================================================

compute_statistics <- function() {
  for(dir in c("results","plot_data","figures","docs")) dir.create(file.path(root,dir),showWarnings=FALSE)
  s <- rd("data/sentinel_selection.tsv")
  h <- rd("data/sentinel_high_groups.tsv")
  p <- rd("data/distance_matched_pairs.tsv.gz")
  base <- s[s$eligible==1,c("cru_id","chrom","start","end","sentinel_acr")]
  key <- function(cru,type,role,id) paste(cru,type,role,id,sep="|")
  tables <- list(); pair_tables <- list()
  for(contrast in c("sentinel_low","sentinel_high")) {
    high <- contrast=="sentinel_high"
    e <- rd(if(high) "data/feature_effects_high.tsv.gz" else "data/feature_effects.tsv.gz")
    ek <- with(e,key(cru_id,feature_type,role,feature_id))
    stopifnot(!anyDuplicated(ek),all(is.finite(e$delta_log2)))
    q <- p
    q$inside_delta_log2 <- e$delta_log2[match(with(q,key(cru_id,feature_type,"inside",observed_feature)),ek)]
    q$control_delta_log2 <- e$delta_log2[match(with(q,key(cru_id,feature_type,"control",control_feature)),ek)]
    q$effects_available <- is.finite(q$inside_delta_log2)&is.finite(q$control_delta_log2)
    q$contrast <- contrast
    pair_tables[[contrast]] <- q
    q <- q[q$effects_available, ]
    tab <- base; tab$contrast <- contrast
    if(high) {
      hs <- h[match(tab$cru_id,h$cru_id), ]
      tab$n_focal_accessions <- hs$n_high; tab$n_normal_accessions <- hs$n_normal
      tab$focal_threshold <- hs$high_threshold; tab$normal_threshold <- hs$normal_threshold
      tab$groups_valid <- hs$high_groups_valid==1
    } else {
      ls <- s[match(tab$cru_id,s$cru_id), ]
      tab$n_focal_accessions <- ls$n_low; tab$n_normal_accessions <- ls$n_normal
      tab$focal_threshold <- ls$low_threshold; tab$normal_threshold <- ls$normal_threshold
      tab$groups_valid <- TRUE
    }
    for(kind in c("ACR","gene")) {
      qq <- q[q$feature_type==kind, ]
      count <- table(qq$cru_id)
      tab[[paste0("n_pairs_",kind)]] <- as.integer(count[tab$cru_id])
      tab[[paste0("n_pairs_",kind)]][is.na(tab[[paste0("n_pairs_",kind)]])] <- 0L
      for(role in c("inside","control")) {
        med <- vapply(split(qq[[paste0(role,"_delta_log2")]],qq$cru_id),median,numeric(1))
        tab[[paste0(role,"_",kind,"_median_delta")]] <- unname(med[tab$cru_id])
      }
    }
    for(role in c("inside","control")) {
      a <- tab[[paste0(role,"_ACR_median_delta")]]
      g <- tab[[paste0(role,"_gene_median_delta")]]
      available <- tab$groups_valid & is.finite(a)&is.finite(g)
      hit <- if(high) a>0 & g>0 else a<0 & g<0
      tab[[paste0(role,"_synchronous")]] <- ifelse(available,as.integer(hit),NA_integer_)
    }
    tab$complete_paired_controls <- is.finite(tab$inside_synchronous)&is.finite(tab$control_synchronous)
    tables[[contrast]] <- tab
  }
  ids <- Reduce(intersect,lapply(tables,function(t) t$cru_id[t$complete_paired_controls]))
  if(!length(ids)) stop("No common CRUs with complete data in both contrasts")
  eligibility <- base
  for(contrast in names(tables)) {
    tab <- tables[[contrast]]
    tab$included_in_both_contrasts <- tab$cru_id %in% ids
    tables[[contrast]] <- tab
    eligibility[[paste0(contrast,"_complete")]] <- tab$complete_paired_controls
  }
  eligibility$included_in_both_contrasts <- eligibility$cru_id %in% ids
  wr(eligibility,"results/CRU_eligibility_audit.tsv")
  wr(do.call(rbind,tables),"results/per_CRU_contrast_results.tsv")
  wr(do.call(rbind,pair_tables),"results/matched_pair_contrast_effects.tsv.gz")

  B <- 10000L; contrast_seed <- 20261021L; pooled_seed <- 20261022L
  # Refresh legacy parameter records without changing the original random seeds.
  parameters <- rd("docs/analysis_parameters.tsv")
  actual <- c(permutation_seed=contrast_seed,pooled_permutation_seed=pooled_seed,permutations=B)
  parameters <- parameters[!parameters$parameter %in% names(actual), ]
  wr(rbind(parameters,data.frame(parameter=names(actual),value=unname(actual))),"docs/analysis_parameters.tsv")
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(contrast_seed)
  tests <- list(); nulls <- list(); source <- list()
  for(contrast in names(tables)) {
    tab <- tables[[contrast]]; tab <- tab[match(ids,tab$cru_id), ]
    x <- tab$inside_synchronous; y <- tab$control_synchronous; d <- x-y
    K <- sum(d!=0); positive <- sum(d==1); observed <- sum(d)
    null_sum <- 2*rbinom(B,K,.5)-K
    log_one <- if(K) pbinom(positive-1,K,.5,lower.tail=FALSE,log.p=TRUE) else 0
    log_two <- if(K) min(0,log(2)+min(log_one,pbinom(positive,K,.5,log.p=TRUE))) else 0
    one <- exp(log_one); two <- exp(log_two)
    if(one==0) one <- NA_real_
    if(two==0) two <- NA_real_
    tests[[contrast]] <- data.frame(contrast=contrast,n_CRUs=length(ids),
      inside_count=sum(x),control_count=sum(y),inside_fraction=mean(x),control_fraction=mean(y),
      difference=mean(d),inside_only=sum(d==1),control_only=sum(d== -1),
      both=sum(x==1 & y==1),neither=sum(x==0 & y==0),permutations=B,
      permutation_P_two_sided=(1+sum(abs(null_sum)>=abs(observed)))/(B+1),
      permutation_P_inside_greater=(1+sum(null_sum>=observed))/(B+1),
      exact_P_two_sided=two,exact_log10_P_two_sided=log_two/log(10),
      exact_P_inside_greater=one,exact_log10_P_inside_greater=log_one/log(10))
    nulls[[contrast]] <- data.frame(contrast=contrast,iteration=seq_len(B),null_sum=null_sum)
    source[[contrast]] <- data.frame(contrast=contrast,
      endpoint=if(contrast=="sentinel_low") "synchronous_decrease" else "synchronous_increase",
      group=c("Inside CRU","Distance-matched controls"),n_CRUs=length(ids),
      n_synchronous=c(sum(x),sum(y)),n_other=c(length(ids)-sum(x),length(ids)-sum(y)),
      synchronous_fraction=c(mean(x),mean(y)),other_fraction=c(1-mean(x),1-mean(y)))
  }
  tests <- do.call(rbind,tests); rownames(tests) <- NULL

  # Fig3 pools the TWO contrast observations per CRU, not unique-CRU outcomes.
  # Keep their dependence: flip both contrast differences together within CRU.
  low <- tables$sentinel_low[match(ids,tables$sentinel_low$cru_id), ]
  high <- tables$sentinel_high[match(ids,tables$sentinel_high$cru_id), ]
  score <- data.frame(cru_id=ids,
    inside_decrease=low$inside_synchronous,inside_increase=high$inside_synchronous,
    control_decrease=low$control_synchronous,control_increase=high$control_synchronous)
  score$inside_sum <- score$inside_decrease+score$inside_increase
  score$control_sum <- score$control_decrease+score$control_increase
  score$difference_sum <- score$inside_sum-score$control_sum
  wr(score,"results/pooled_per_CRU_scores.tsv")
  n1 <- sum(abs(score$difference_sum)==1); n2 <- sum(abs(score$difference_sum)==2)
  observed <- sum(score$difference_sum)
  set.seed(pooled_seed)
  pooled_null <- (2*rbinom(B,n1,.5)-n1)+2*(2*rbinom(B,n2,.5)-n2)
  # Exact CRU-clustered null: sum of n1 unit-weight and n2 double-weight signs.
  # Integrate over Binomial(n2,.5), using log probabilities to prevent underflow.
  weighted_tail_log <- function(threshold,n1,n2) {
    j <- 0:n2
    required <- ceiling((threshold-(4*j-2*n2)+n1)/2)
    terms <- dbinom(j,n2,.5,log=TRUE)+pbinom(required-1,n1,.5,lower.tail=FALSE,log.p=TRUE)
    mx <- max(terms)
    if(!is.finite(mx)) return(-Inf)
    mx+log(sum(exp(terms-mx)))
  }
  pooled_log_one <- weighted_tail_log(observed,n1,n2)
  pooled_log_two <- if(observed==0) 0 else min(0,log(2)+weighted_tail_log(abs(observed),n1,n2))
  pooled_exact <- exp(pooled_log_two); if(pooled_exact==0) pooled_exact <- NA_real_
  pooled_one <- exp(pooled_log_one); if(pooled_one==0) pooled_one <- NA_real_
  pooled <- data.frame(endpoint="synchronous_change_pooled",n_unique_CRUs=length(ids),
    contrasts_per_CRU=2L,n_observations=2L*length(ids),inside_count=sum(score$inside_sum),
    control_count=sum(score$control_sum),inside_other=sum(2L-score$inside_sum),
    control_other=sum(2L-score$control_sum),inside_fraction=mean(score$inside_sum)/2,
    control_fraction=mean(score$control_sum)/2,difference=mean(score$difference_sum)/2,
    n_CRUs_abs_difference_1=n1,n_CRUs_abs_difference_2=n2,observed_sum=observed,
    permutations=B,permutation_unit="CRU_both_contrasts_together",
    permutation_P_two_sided=(1+sum(abs(pooled_null)>=abs(observed)))/(B+1),
    permutation_P_inside_greater=(1+sum(pooled_null>=observed))/(B+1),
    exact_P_two_sided=pooled_exact,exact_log10_P_two_sided=pooled_log_two/log(10),
    exact_P_inside_greater=pooled_one,exact_log10_P_inside_greater=pooled_log_one/log(10))
  family_mc <- p.adjust(c(tests$permutation_P_two_sided,pooled$permutation_P_two_sided),"holm")
  family_exact <- p.adjust(c(tests$exact_P_two_sided,pooled$exact_P_two_sided),"holm")
  tests$permutation_Holm_P_two_sided <- family_mc[1:2]
  tests$exact_Holm_P_two_sided <- family_exact[1:2]
  pooled$permutation_Holm_P_two_sided <- family_mc[3]
  pooled$exact_Holm_P_two_sided <- family_exact[3]
  wr(pooled,"results/pooled_synchronous_change_test.tsv")
  wr(data.frame(iteration=seq_len(B),null_sum=pooled_null),"results/pooled_CRU_sign_flip_null.tsv.gz")
  pooled_source <- do.call(rbind,lapply(c("Inside CRU","Distance-matched controls"),function(group) {
    internal <- group=="Inside CRU"
    decrease <- sum(if(internal) score$inside_decrease else score$control_decrease)
    increase <- sum(if(internal) score$inside_increase else score$control_increase)
    n <- 2L*length(ids)
    data.frame(endpoint="synchronous_change_pooled",group=group,n_unique_CRUs=length(ids),
      n_observations=n,n_decrease=decrease,n_increase=increase,
      n_not_decreased=length(ids)-decrease,n_not_increased=length(ids)-increase,
      n_synchronous=decrease+increase,n_other=n-decrease-increase,
      synchronous_fraction=(decrease+increase)/n,other_fraction=(n-decrease-increase)/n,
      permutation_P_two_sided=pooled$permutation_P_two_sided,
      permutation_Holm_P_two_sided=pooled$permutation_Holm_P_two_sided,
      exact_P_two_sided=pooled$exact_P_two_sided,exact_Holm_P_two_sided=pooled$exact_Holm_P_two_sided)
  }))
  wr(pooled_source,"plot_data/pooled_synchronous_change.tsv")
  wr(tests,"results/paired_contrast_tests.tsv")
  wr(do.call(rbind,nulls),"results/contrast_permutation_null.tsv.gz")
  source <- do.call(rbind,source); rownames(source) <- NULL
  idx <- match(source$contrast,tests$contrast)
  for(col in c("permutation_P_two_sided","permutation_Holm_P_two_sided","exact_P_two_sided","exact_Holm_P_two_sided")) {
    source[[col]] <- tests[[col]][idx]
  }
  wr(source,"plot_data/sentinel_contrast_fractions.tsv")
  summary <- c("SEPARATE SENTINEL-LOW AND SENTINEL-HIGH CONTRASTS",
    "Same original sentinel, same frozen distance-matched feature pairs, equal pair counts inside/control.",
    "Low: sentinel <=Q20 vs >=Q50; synchronous decrease if both feature-type medians <0.",
    "High: sentinel >=Q80 vs <=Q50; synchronous increase if both feature-type medians >0.",
    "Effects are mean(log2(raw+1),focal)-mean(log2(raw+1),normal).",
    "Only CRUs with complete paired data in BOTH contrasts enter either figure.",
    "Fig1/2 use paired sign flips; Fig3 pools counts and flips the TWO observations together within CRU.",
    sprintf("Holm correction covers THREE panels; Fig3 has %d observations from %d CRUs (two contrasts per CRU).",2L*length(ids),length(ids)),
    "")
  for(i in seq_len(nrow(tests))) summary <- c(summary,sprintf(
    "%s: inside %d/%d (%.4f%%); control %d/%d (%.4f%%); two-sided permutation P=%.8g; Holm P=%.8g; exact two-sided P=%.8g",
    tests$contrast[i],tests$inside_count[i],tests$n_CRUs[i],100*tests$inside_fraction[i],
    tests$control_count[i],tests$n_CRUs[i],100*tests$control_fraction[i],
    tests$permutation_P_two_sided[i],tests$permutation_Holm_P_two_sided[i],tests$exact_P_two_sided[i]))
  summary <- c(summary,sprintf(
    "Fig3 synchronous change (pooled counts): inside %d/%d (%.4f%%), other %d; controls %d/%d (%.4f%%), other %d; CRU-clustered two-sided permutation P=%.8g; Holm P=%.8g; exact clustered P=%.8g",
    pooled$inside_count,pooled$n_observations,100*pooled$inside_fraction,pooled$inside_other,
    pooled$control_count,pooled$n_observations,100*pooled$control_fraction,pooled$control_other,
    pooled$permutation_P_two_sided,pooled$permutation_Holm_P_two_sided,pooled$exact_P_two_sided),
    "Other means NOT synchronous for that contrast, not evidence of zero accessibility/expression change.")
  writeLines(summary,file.path(root,"docs/RESULTS_SUMMARY.txt"))
  capture.output(sessionInfo(),file=file.path(root,"docs/R_sessionInfo.txt"))
  cat(paste(summary,collapse="\n"),"\n")
}


# ========================================================================
# plot_figures
# ========================================================================

plot_figures <- function() {
  figdir <- file.path(root,"figures")
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  d <- read.delim(file.path(root,"plot_data/sentinel_contrast_fractions.tsv"),stringsAsFactors=FALSE)
  FONT_PT <- 7; LINE_PT <- .75; LWD <- LINE_PT/(72/96)
  FONT_FAMILY <- Sys.getenv("CRU_FONT_FAMILY", unset="Arial")
  open_pdf <- function(filename) {
    if(capabilities("cairo")) {
      cairo_pdf(filename,width=5.1,height=4.1,family=FONT_FAMILY,pointsize=FONT_PT)
    } else {
      pdf(filename,width=5.1,height=4.1,family="Helvetica",pointsize=FONT_PT)
    }
  }
  if(!capabilities("cairo")) FONT_FAMILY <- "Helvetica"
  draw <- function(contrast,filename,title,color) {
    x <- d[d$contrast==contrast, ]; x <- x[match(c("Inside CRU","Distance-matched controls"),x$group), ]
    stopifnot(nrow(x)==2L,all(is.finite(x$synchronous_fraction)))
    decreasing <- contrast=="sentinel_low"
    verb <- if(decreasing) "decrease" else "increase"
    open_pdf(file.path(figdir,filename))
    on.exit(dev.off())
    par(family=FONT_FAMILY,cex=1,las=1,lwd=LWD,fg="#202020",col="#202020",
        mar=c(6.3,5.2,3.8,1.0),mgp=c(2.5,.65,0))
    values <- rbind(x$synchronous_fraction,x$other_fraction)
    bp <- barplot(values,col=c(color,"#D8D8D8"),border="white",lwd=LWD,
                  names.arg=rep("",2),ylim=c(0,1),axes=FALSE)
    axis(2,at=seq(0,1,.2),labels=sprintf("%.1f",seq(0,1,.2)),lwd=LWD,lwd.ticks=LWD,tck=-.025,cex.axis=1)
    mtext(x$group,side=1,at=bp,line=1.2,cex=1)
    title(ylab="Fraction of CRUs",line=2.8,cex.lab=1)
    mtext(title,side=3,line=2.2,font=2,cex=1)
    mtext(sprintf("Paired permutation (two-sided): P = %.2g; Holm P = %.2g",
      x$permutation_P_two_sided[1],x$permutation_Holm_P_two_sided[1]),side=3,line=.8,cex=1)
    for(i in 1:2) {
      bottom <- if(values[1,i]>=.17) paste0("synchronous\n",verb,"\n",sprintf("%.1f%%",100*values[1,i])) else
        sprintf("%.1f%%",100*values[1,i])
      top <- if(values[2,i]>=.17) paste0("not both\n",if(decreasing) "decreased" else "increased","\n",
                                      sprintf("%.1f%%",100*values[2,i])) else sprintf("%.1f%%",100*values[2,i])
      text(bp[i],values[1,i]/2,bottom,col="white",font=2,cex=1)
      text(bp[i],values[1,i]+values[2,i]/2,top,col="#4D4D4D",cex=1)
    }
    mtext(if(decreasing) "Sentinel lowest 20% vs upper 50%" else "Sentinel highest 20% vs lower 50%",
          side=1,line=3.3,cex=1)
    mtext(paste0("Equal-sized matched features; same complete paired CRUs: n=",x$n_CRUs[1]),
          side=1,line=4.7,cex=1)
  }
  draw("sentinel_low","Fig1_CRU_synchronous_decrease.pdf",
       "Sentinel-low: other ACRs and genes decrease together","#3F8B5B")
  draw("sentinel_high","Fig2_CRU_synchronous_increase.pdf",
       "Sentinel-high: other ACRs and genes increase together","#B55C31")

  draw_pooled <- function() {
    x <- read.delim(file.path(root,"plot_data/pooled_synchronous_change.tsv"),stringsAsFactors=FALSE)
    x <- x[match(c("Inside CRU","Distance-matched controls"),x$group), ]
    stopifnot(nrow(x)==2L,all(is.finite(x$synchronous_fraction)))
    open_pdf(file.path(figdir,"Fig3_CRU_synchronous_change.pdf"))
    on.exit(dev.off())
    par(family=FONT_FAMILY,cex=1,las=1,lwd=LWD,fg="#202020",col="#202020",
        mar=c(6.3,5.2,3.8,1.0),mgp=c(2.5,.65,0))
    values <- rbind(x$synchronous_fraction,x$other_fraction)
    bp <- barplot(values,col=c("#3F8B5B","#D8D8D8"),border="white",lwd=LWD,
                  names.arg=rep("",2),ylim=c(0,1),axes=FALSE)
    axis(2,at=seq(0,1,.2),labels=sprintf("%.1f",seq(0,1,.2)),lwd=LWD,lwd.ticks=LWD,tck=-.025,cex.axis=1)
    mtext(x$group,side=1,at=bp,line=1.2,cex=1)
    title(ylab="Fraction of CRU-contrast observations",line=2.8,cex.lab=1)
    mtext("Fig3  Synchronous changes (decrease + increase)",side=3,line=2.2,font=2,cex=1)
    mtext(sprintf("CRU-clustered sign-flip: P = %.2g; Holm P = %.2g",
      x$permutation_P_two_sided[1],x$permutation_Holm_P_two_sided[1]),side=3,line=.8,cex=1)
    text(bp,values[1,]/2,paste0("synchronous\nchange\n",sprintf("%.1f%%",100*values[1,])),
         col="white",font=2,cex=1)
    text(bp,values[1,]+values[2,]/2,paste0("not synchronous\n",sprintf("%.1f%%",100*values[2,])),
         col="#4D4D4D",cex=1)
    mtext("Decrease + increase counts pooled; each CRU contributes twice",side=1,line=3.3,cex=1)
    mtext(sprintf("%s observations (%s CRUs x 2 contrasts)",
          format(x$n_observations[1],big.mark=",",trim=TRUE),format(x$n_unique_CRUs[1],big.mark=",",trim=TRUE)),
          side=1,line=4.7,cex=1)
  }
  draw_pooled()
  cat("DONE: separate low/high PDFs and Fig3 pooled synchronous-change PDF.\n")
}


# ========================================================================
# validate_results
# ========================================================================

validate_results <- function() {
  checks <- list()
  check <- function(ok,label) {
    checks[[length(checks)+1L]] <<- data.frame(check=label,status=if(isTRUE(ok)) "PASS" else "FAIL")
    if(!isTRUE(ok)) stop("Validation failed: ",label)
  }
  near <- function(x,y,label,tol=1e-10) check(length(x)==length(y) && all(is.finite(x)) &&
    all(is.finite(y)) && max(abs(x-y),0)<tol,label)
  s <- rd("data/sentinel_selection.tsv"); h <- rd("data/sentinel_high_groups.tsv")
  z <- rd("data/sentinel_accessibility.tsv.gz"); a <- rd("data/accessions.tsv")
  p <- rd("data/distance_matched_pairs.tsv.gz")
  t <- rd("results/per_CRU_contrast_results.tsv"); tests <- rd("results/paired_contrast_tests.tsv")
  nul <- rd("results/contrast_permutation_null.tsv.gz")
  source <- rd("plot_data/sentinel_contrast_fractions.tsv")
  audit <- rd("results/CRU_eligibility_audit.tsv")
  parameters <- rd("docs/analysis_parameters.tsv")
  param <- function(k) as.numeric(parameters$value[match(k,parameters$parameter)])
  cru_path <- file.path(root,"data/cru_input.bed")
  raw_cru <- read.delim(cru_path,header=FALSE,stringsAsFactors=FALSE)
  check(!anyDuplicated(raw_cru[[4]]) && setequal(s$cru_id,raw_cru[[4]]),"new_input_CRU_identifiers_preserved")
  cru_idx <- match(s$cru_id,raw_cru[[4]])
  check(all(s$chrom==sub("^Chr","",raw_cru[[1]][cru_idx])) &&
        all(s$start==raw_cru[[2]][cru_idx]) && all(s$end==raw_cru[[3]][cru_idx]),"new_input_CRU_coordinates_preserved")
  check(nrow(s)==param("n_input_CRUs") && sum(s$has_sentinel)==param("n_sentinels") &&
        sum(s$eligible)==param("n_eligible_CRUs"),"new_CRUs_and_sentinels_match_preparation")
  check(nrow(a)==200L && identical(names(z)[-1],a$accession),"same_200_accessions_and_alignment")
  check(!anyDuplicated(s$cru_id) && !anyDuplicated(h$cru_id) && !anyDuplicated(z$cru_id),"unique_sentinel_identifiers")
  check(nrow(z)==sum(s$has_sentinel) && setequal(z$cru_id,h$cru_id) &&
        setequal(z$cru_id,s$cru_id[s$has_sentinel==1]),"all_sentinel_raw_values_and_high_groups_recorded")
  mat <- as.matrix(z[,-1]); hs <- h[match(z$cru_id,h$cru_id), ]; ls <- s[match(z$cru_id,s$cru_id), ]
  near(apply(mat,1,quantile,probs=.8,type=7,na.rm=TRUE),hs$high_threshold,"high_threshold_is_raw_Q80")
  near(apply(mat,1,quantile,probs=.5,type=7,na.rm=TRUE),hs$normal_threshold,"high_normal_threshold_is_raw_Q50")
  near(hs$normal_threshold,ls$normal_threshold,"same_sentinel_median_boundary")
  hi <- mat>=hs$high_threshold; normal_hi <- mat<=hs$normal_threshold
  lo <- mat<=ls$low_threshold; normal_lo <- mat>=ls$normal_threshold
  near(rowSums(hi,na.rm=TRUE),hs$n_high,"high_group_accession_counts")
  near(rowSums(normal_hi,na.rm=TRUE),hs$n_normal,"lower_half_normal_accession_counts")
  near(rowSums(lo,na.rm=TRUE),ls$n_low,"low_group_counts_unchanged")
  near(rowSums(normal_lo,na.rm=TRUE),ls$n_normal,"upper_half_normal_counts_unchanged")
  valid_hi <- hs$high_groups_valid==1
  check(all(hs$high_threshold[valid_hi]>hs$normal_threshold[valid_hi]) &&
        all(rowSums(hi & normal_hi,na.rm=TRUE)[valid_hi]==0),"high_and_normal_groups_do_not_overlap")
  check(all(hs$n_high[valid_hi]>=15 & hs$n_normal[valid_hi]>=15),"high_group_minimum_size")
  means <- function(mask) rowSums(log2(mat+1)*mask,na.rm=TRUE)/rowSums(mask,na.rm=TRUE)
  near((means(hi)-means(normal_hi))[valid_hi],hs$sentinel_delta_high_vs_normal[valid_hi],"sentinel_high_delta_from_raw_values")
  check(all(hs$sentinel_delta_high_vs_normal[valid_hi]>0),"sentinel_accessibility_really_increases")
  check(sum(p$feature_type=="ACR")==param("n_ACR_pairs") && sum(p$feature_type=="gene")==param("n_gene_pairs"),"genomic_matches_match_preparation")
  check(all(p$control_cru_id!=p$cru_id),"controls_outside_focal_CRU")
  near(p$tolerance_bp,pmax(10000,floor(.2*p$observed_distance_bp)),"matching_tolerance")
  check(all(p$absolute_distance_error_bp<=p$tolerance_bp),"control_distances_within_tolerance")
  key <- function(cru,type,role,id) paste(cru,type,role,id,sep="|")
  for(contrast in c("sentinel_low","sentinel_high")) {
    high <- contrast=="sentinel_high"
    e <- rd(if(high) "data/feature_effects_high.tsv.gz" else "data/feature_effects.tsv.gz")
    focal <- if(high) "mean_high_log2" else "mean_low_log2"
    near(e$delta_log2,e[[focal]]-e$mean_normal_log2,paste0(contrast,"_feature_effects_from_group_means"))
    check(all(e[[if(high) "n_high" else "n_low"]]>=15 & e$n_normal>=15),paste0(contrast,"_feature_sample_sizes"))
    ek <- with(e,key(cru_id,feature_type,role,feature_id))
    check(!anyDuplicated(ek),paste0(contrast,"_unique_feature_effects"))
    check(!any(e$feature_type=="ACR" & e$role=="inside" &
              e$feature_id==s$sentinel_acr[match(e$cru_id,s$cru_id)]),paste0(contrast,"_sentinel_excluded"))
    tab <- t[t$contrast==contrast, ]
    check(nrow(tab)==sum(s$eligible) && !anyDuplicated(tab$cru_id) &&
          setequal(tab$cru_id,s$cru_id[s$eligible==1]),paste0(contrast,"_eligible_CRUs_accounted_for"))
    grouping <- if(high) h[match(tab$cru_id,h$cru_id), ] else s[match(tab$cru_id,s$cru_id), ]
    near(tab$n_focal_accessions,grouping[[if(high) "n_high" else "n_low"]],paste0(contrast,"_recorded_focal_accession_counts"))
    near(tab$n_normal_accessions,grouping$n_normal,paste0(contrast,"_recorded_normal_accession_counts"))
    near(tab$focal_threshold,grouping[[if(high) "high_threshold" else "low_threshold"]],paste0(contrast,"_recorded_focal_thresholds"))
    near(tab$normal_threshold,grouping$normal_threshold,paste0(contrast,"_recorded_normal_thresholds"))
    q <- p
    q$inside <- e$delta_log2[match(with(q,key(cru_id,feature_type,"inside",observed_feature)),ek)]
    q$control <- e$delta_log2[match(with(q,key(cru_id,feature_type,"control",control_feature)),ek)]
    q <- q[is.finite(q$inside)&is.finite(q$control), ]
    for(kind in c("ACR","gene")) {
      qq <- q[q$feature_type==kind, ]; count <- table(qq$cru_id)
      nc <- as.integer(count[tab$cru_id]); nc[is.na(nc)] <- 0L
      near(tab[[paste0("n_pairs_",kind)]],nc,paste0(contrast,"_",kind,"_equal_pair_counts"))
      for(role in c("inside","control")) {
        med <- vapply(split(qq[[role]],qq$cru_id),median,numeric(1))
        idx <- match(tab$cru_id,names(med)); have <- !is.na(idx)
        near(tab[[paste0(role,"_",kind,"_median_delta")]][have],unname(med[idx[have]]),
             paste(contrast,role,kind,"pair_medians",sep="_"))
        check(all(is.na(tab[[paste0(role,"_",kind,"_median_delta")]][!have])),
              paste(contrast,role,kind,"missing_medians_stay_NA",sep="_"))
      }
    }
    for(role in c("inside","control")) {
      aa <- tab[[paste0(role,"_ACR_median_delta")]]; gg <- tab[[paste0(role,"_gene_median_delta")]]
      have <- tab$groups_valid & is.finite(aa)&is.finite(gg)
      hit <- if(high) aa[have]>0 & gg[have]>0 else aa[have]<0 & gg[have]<0
      near(tab[[paste0(role,"_synchronous")]][have],as.integer(hit),paste(contrast,role,"synchronous_definition",sep="_"))
      check(all(is.na(tab[[paste0(role,"_synchronous")]][!have])),paste(contrast,role,"missing_is_not_zero",sep="_"))
    }
    available <- is.finite(tab$inside_synchronous)&is.finite(tab$control_synchronous)
    check(identical(tab$complete_paired_controls,available),paste0(contrast,"_complete_pair_flag"))
    included <- tab$included_in_both_contrasts
    check(all(tab$complete_paired_controls[included]),paste0(contrast,"_only_complete_CRUs_in_figures"))
    check(setequal(tab$cru_id[included],audit$cru_id[audit$included_in_both_contrasts]),paste0(contrast,"_common_CRU_set"))
    x <- tab$inside_synchronous[included]; y <- tab$control_synchronous[included]; dd <- x-y
    tt <- tests[tests$contrast==contrast, ]; nu <- nul[nul$contrast==contrast, ]; fig <- source[source$contrast==contrast, ]
    check(nrow(tt)==1L && nrow(fig)==2L,paste0(contrast,"_one_test_two_bars"))
    near(c(tt$n_CRUs,tt$inside_count,tt$control_count),c(length(x),sum(x),sum(y)),paste0(contrast,"_test_counts"))
    near(c(tt$inside_fraction,tt$control_fraction,tt$difference),c(mean(x),mean(y),mean(dd)),paste0(contrast,"_test_fractions"))
    near(c(tt$inside_only,tt$control_only,tt$both,tt$neither),
         c(sum(dd==1),sum(dd== -1),sum(x==1 & y==1),sum(x==0 & y==0)),paste0(contrast,"_paired_contingency_cells"))
    check(nrow(nu)==10000L && identical(nu$iteration,seq_len(10000L)),paste0(contrast,"_10000_saved_null_draws"))
    near(tt$permutation_P_two_sided,(1+sum(abs(nu$null_sum)>=abs(sum(dd))))/10001,paste0(contrast,"_two_sided_permutation_P"),1e-12)
    near(tt$permutation_P_inside_greater,(1+sum(nu$null_sum>=sum(dd)))/10001,paste0(contrast,"_one_sided_permutation_P"),1e-12)
    K <- sum(dd!=0)
    exact <- if(K) binom.test(sum(dd==1),K,.5,alternative="two.sided")$p.value else 1
    near(tt$exact_P_two_sided,exact,paste0(contrast,"_independent_exact_binom_test"),1e-12)
    fig <- fig[match(c("Inside CRU","Distance-matched controls"),fig$group), ]
    near(fig$n_synchronous,c(sum(x),sum(y)),paste0(contrast,"_plot_counts"))
    near(fig$synchronous_fraction,c(mean(x),mean(y)),paste0(contrast,"_plot_fractions"))
    near(fig$synchronous_fraction+fig$other_fraction,c(1,1),paste0(contrast,"_stacked_bars_sum_to_one"))
    near(fig$n_synchronous+fig$n_other,rep(length(x),2),paste0(contrast,"_plot_denominators"))
    near(fig$permutation_P_two_sided,rep(tt$permutation_P_two_sided,2),paste0(contrast,"_plot_P_values"))
    near(fig$permutation_Holm_P_two_sided,rep(tt$permutation_Holm_P_two_sided,2),paste0(contrast,"_plot_Holm_P_values"))
  }
  check(sum(audit$included_in_both_contrasts)>0L &&
        all(tests$n_CRUs==sum(audit$included_in_both_contrasts)),"same_complete_CRUs_in_both_figures")
  check(identical(audit$included_in_both_contrasts,audit$sentinel_low_complete & audit$sentinel_high_complete),
        "shared_subset_is_exact_intersection_of_complete_contrasts")
  low_tab <- t[t$contrast=="sentinel_low" & t$included_in_both_contrasts, ]
  high_tab <- t[t$contrast=="sentinel_high" & t$included_in_both_contrasts, ]
  high_tab <- high_tab[match(low_tab$cru_id,high_tab$cru_id), ]
  check(identical(low_tab$sentinel_acr,high_tab$sentinel_acr),"same_sentinel_in_low_and_high_figures")
  near(low_tab$n_pairs_ACR,high_tab$n_pairs_ACR,"same_ACR_pair_counts_across_contrasts")
  near(low_tab$n_pairs_gene,high_tab$n_pairs_gene,"same_gene_pair_counts_across_contrasts")
  lo_test <- tests[tests$contrast=="sentinel_low", ]
  check(lo_test$inside_count==sum(low_tab$inside_synchronous) &&
        lo_test$control_count==sum(low_tab$control_synchronous),"strict_matched_decrease_counts_from_new_CRUs")
  pooled <- rd("results/pooled_synchronous_change_test.tsv")
  pooled_source <- rd("plot_data/pooled_synchronous_change.tsv")
  score <- rd("results/pooled_per_CRU_scores.tsv")
  pooled_null <- rd("results/pooled_CRU_sign_flip_null.tsv.gz")
  check(nrow(pooled)==1L && nrow(pooled_source)==2L,"one_pooled_test_two_pooled_bars")
  check(identical(score$cru_id,low_tab$cru_id),"pooled_scores_preserve_common_CRU_alignment")
  near(score$inside_decrease,low_tab$inside_synchronous,"pooled_inside_decrease_from_Fig1")
  near(score$inside_increase,high_tab$inside_synchronous,"pooled_inside_increase_from_Fig2")
  near(score$control_decrease,low_tab$control_synchronous,"pooled_control_decrease_from_Fig1")
  near(score$control_increase,high_tab$control_synchronous,"pooled_control_increase_from_Fig2")
  near(score$inside_sum,score$inside_decrease+score$inside_increase,"pooled_inside_counts_sum_directions")
  near(score$control_sum,score$control_decrease+score$control_increase,"pooled_control_counts_sum_directions")
  near(score$difference_sum,score$inside_sum-score$control_sum,"pooled_paired_difference_per_CRU")
  near(pooled$n_unique_CRUs,nrow(score),"pooled_unique_CRUs")
  near(pooled$n_observations,2*nrow(score),"pooled_denominator_is_twice_CRUs")
  check(pooled$n_observations==2L*nrow(low_tab) && pooled$n_unique_CRUs==nrow(low_tab),"pooled_observations_not_independent_CRUs")
  near(c(pooled$inside_count,pooled$control_count,pooled$inside_other,pooled$control_other),
       c(sum(score$inside_sum),sum(score$control_sum),sum(2-score$inside_sum),sum(2-score$control_sum)),"pooled_synchronous_and_other_counts")
  near(c(pooled$inside_count,pooled$control_count,pooled$inside_other,pooled$control_other),
       c(sum(tests$inside_count),sum(tests$control_count),
         sum(tests$n_CRUs-tests$inside_count),sum(tests$n_CRUs-tests$control_count)),"Fig3_correct_pooled_counts")
  near(c(pooled$inside_fraction,pooled$control_fraction),
       c(mean(score$inside_sum)/2,mean(score$control_sum)/2),"pooled_fractions_use_pooled_denominator")
  ps <- pooled_source[match(c("Inside CRU","Distance-matched controls"),pooled_source$group), ]
  near(ps$n_synchronous,c(pooled$inside_count,pooled$control_count),"Fig3_source_synchronous_counts")
  near(ps$n_other,c(pooled$inside_other,pooled$control_other),"Fig3_source_other_counts")
  near(ps$n_synchronous,ps$n_decrease+ps$n_increase,"Fig3_colored_segment_sums_direction_counts")
  near(ps$n_other,ps$n_not_decreased+ps$n_not_increased,"Fig3_gray_segment_sums_non_synchronous_counts")
  near(ps$n_synchronous+ps$n_other,ps$n_observations,"Fig3_segments_share_pooled_denominator")
  near(ps$synchronous_fraction,c(pooled$inside_fraction,pooled$control_fraction),"Fig3_source_fractions")
  near(ps$synchronous_fraction+ps$other_fraction,c(1,1),"Fig3_bars_sum_to_one")
  near(ps$synchronous_fraction,(tests$inside_fraction[1]+tests$inside_fraction[2])/2*c(1,0)+
       (tests$control_fraction[1]+tests$control_fraction[2])/2*c(0,1),"Fig3_fraction_is_average_not_sum_of_percentages")
  check(pooled$permutation_unit=="CRU_both_contrasts_together","Fig3_permutation_preserves_within_CRU_dependence")
  n1 <- sum(abs(score$difference_sum)==1); n2 <- sum(abs(score$difference_sum)==2)
  near(c(pooled$n_CRUs_abs_difference_1,pooled$n_CRUs_abs_difference_2),c(n1,n2),"pooled_cluster_sign_weights")
  check(nrow(pooled_null)==10000L && identical(pooled_null$iteration,seq_len(10000L)),"Fig3_10000_saved_clustered_draws")
  check(all(abs(pooled_null$null_sum)<=n1+2*n2) && all((pooled_null$null_sum-n1)%%2==0),"Fig3_clustered_null_support_and_parity")
  near(pooled$observed_sum,sum(score$difference_sum),"Fig3_observed_clustered_sum")
  near(pooled$permutation_P_two_sided,
       (1+sum(abs(pooled_null$null_sum)>=abs(sum(score$difference_sum))))/10001,"Fig3_clustered_two_sided_permutation_P",1e-12)
  # Independent exact check by sequential probability convolution, not the
  # binomial-tail integration used by the analysis script.
  distribution <- 1
  for(weight in c(rep(1L,n1),rep(2L,n2))) {
    ii <- seq_along(distribution)
    next_distribution <- numeric(length(distribution)+2L*weight)
    next_distribution[ii] <- distribution/2
    next_distribution[ii+2L*weight] <- next_distribution[ii+2L*weight]+distribution/2
    distribution <- next_distribution
  }
  support <- seq.int(-(n1+2*n2),n1+2*n2)
  near(sum(distribution),1,"Fig3_exact_clustered_distribution_normalized")
  independent_exact <- sum(distribution[abs(support)>=abs(sum(score$difference_sum))])
  near(log10(independent_exact),pooled$exact_log10_P_two_sided,"Fig3_exact_clustered_P_independent_convolution",1e-9)
  mc_family <- p.adjust(c(tests$permutation_P_two_sided,pooled$permutation_P_two_sided),"holm")
  ex_family <- p.adjust(c(tests$exact_P_two_sided,pooled$exact_P_two_sided),"holm")
  near(c(tests$permutation_Holm_P_two_sided,pooled$permutation_Holm_P_two_sided),mc_family,"Holm_correction_THREE_panels")
  near(c(tests$exact_Holm_P_two_sided,pooled$exact_Holm_P_two_sided),ex_family,"exact_Holm_correction_THREE_panels")
  near(ps$permutation_P_two_sided,rep(pooled$permutation_P_two_sided,2),"Fig3_source_clustered_P")
  near(ps$permutation_Holm_P_two_sided,rep(pooled$permutation_Holm_P_two_sided,2),"Fig3_source_Holm_P")
  check(nrow(tests)==2L && setequal(tests$contrast,c("sentinel_low","sentinel_high")),"only_two_genuine_contrast_tests")
  pdfs <- c("Fig1_CRU_synchronous_decrease.pdf","Fig2_CRU_synchronous_increase.pdf","Fig3_CRU_synchronous_change.pdf")
  check(setequal(list.files(file.path(root,"figures")),pdfs),"only_three_requested_PDF_figures")
  for(file in pdfs) {
    check(file.info(file.path(root,"figures",file))$size>0,paste0(file,"_nonempty"))
  }
  write.table(do.call(rbind,checks),file.path(root,"docs/validation_report.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  cat("PASS:",length(checks),"low/high grouping, raw summaries, matching, paired statistics and figure checks.\n")
}

if(args[1] %in% c("all","statistics")) compute_statistics()
if(args[1] %in% c("all","plot")) plot_figures()
if(args[1] %in% c("all","validate")) validate_results()
