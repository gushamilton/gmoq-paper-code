#!/usr/bin/env Rscript
# Summary tables of simulation studies 1-6
#
# What: reads the outputs of the four runs (studies 1-2, study 3, studies 4-5, study 6) and
#   writes the per-study summaries imported by this repository: Table S1 (design of each
#   study), Table S2 (numerical summary per cell: bias, RMSE, mean SE, empirical SD, coverage,
#   SE ratios, timings), and the means over the three timing repeats of studies 2 and 5. It
#   also copies the per-run CSVs beside them and draws working versions of the figures; the
#   paper's Supplementary Figs. S1-S6 are drawn by figures/figS01_*.R ... figS06_*.R from the
#   imported summaries.
# Paper: Supplementary Tables S1, S2; Supplementary Figs. S1-S6 (studies 1-6)
# Input: <ROOT>/.local/ run directories of the four runs (paths in sections below)
# Output: <ROOT>/outputs/simulation_results/{TableS1_design_provenance,
#   TableS2_numerical_summary,S2_mean_timings,S5_mean_timings}.csv (copied to
#   ../../../inputs/simulation_summaries/), plus copies of the run CSVs and working plots
# Runs: on a cluster, under a minute, from <ROOT>:
#   Rscript --vanilla validation/figures/plot_simulation_supplement.R outputs/simulation_results
#   The run directory names below are those of the authors' runs.
a <- commandArgs(TRUE); if (length(a) != 1L) stop("Usage: plot_simulation_supplement.R OUTPUT_DIR")
sc <- grep("^--file=", commandArgs(), value=TRUE)[1]; root <- normalizePath(file.path(dirname(sub("^--file=", "", sc)), "../.."))
requested <- if (grepl("^/",a[1])) a[1] else file.path(root,a[1])
out <- normalizePath(requested, mustWork=FALSE)
canonical <- normalizePath(file.path(root,"outputs","simulation_results"), mustWork=FALSE)
if (!identical(out, canonical)) stop("renderer accepts only canonical output directory: ", canonical)
generated_assets <- c(paste0(c("S1_beta","S1_SE","S2","S3","S4_beta","S4_SE","S4_truth","S5","S6A","S6B","S6C"),".png"),paste0(c("S1_beta","S1_SE","S2","S3","S4_beta","S4_SE","S4_truth","S5","S6A","S6B","S6C"),".pdf"),c("S1_pairs.csv","S1_SE.csv","S2_timings.csv","S2_mean_timings.csv","S3_results.csv","S4_results.csv","S5_panel_timings.csv","S5_mean_timings.csv","S6_stage1_tau.csv","S6_stage1_summary.csv","S6_stage2_parameters.csv","S6_summary.csv","S6C_effects.csv","S6C_summary.csv","TableS1_design_provenance.csv","TableS2_numerical_summary.csv","captions.txt","README.md","render_metadata.txt"))
host <- Sys.info()[["nodename"]]
if (!requireNamespace("ggplot2", quietly=TRUE) || !requireNamespace("jsonlite", quietly=TRUE)) stop("ggplot2 and jsonlite are required")
suppressPackageStartupMessages(library(ggplot2)); suppressPackageStartupMessages(library(jsonlite))
if (file.exists(out)) {
 if (!dir.exists(out)) stop("canonical output path is not a directory: ", out)
 existing <- list.files(out, all.files=TRUE, no..=TRUE)
 unexpected <- setdiff(existing, generated_assets)
 if (length(unexpected)) stop("refusing canonical directory with unexpected existing assets: ", paste(unexpected, collapse=", "))
 if (length(existing)) unlink(file.path(out, existing), force=TRUE)
} else dir.create(out, recursive=TRUE)
readc <- function(p, cols=character()) { if (!file.exists(p)) stop("missing input: ", p); z <- read.csv(p, stringsAsFactors=FALSE, check.names=FALSE); if (!all(cols %in% names(z))) stop("bad columns in ", p); z }
num <- function(z, cols) { for (x in cols) { z[[x]] <- suppressWarnings(as.numeric(z[[x]])); if (any(!is.finite(z[[x]]))) stop("non-finite ", x) }; z }
savep <- function(p, stem, w=8, h=6) { ggsave(file.path(out,paste0(stem,".png")),p,width=w,height=h,units="in",dpi=300) }
plain <- theme_bw(base_size=10) + theme(legend.position="bottom", panel.grid.minor=element_blank())
rootcsv <- function(p) file.path(root,p)
copycsv <- function(p, n) if (file.exists(p)) file.copy(p,file.path(out,n))
group_rows <- function(d, keys) { key <- do.call(paste,c(lapply(d[keys],as.character),sep="\r")); split(seq_len(nrow(d)),factor(key,levels=unique(key))) }
mean_timings <- function(d, keys, metrics) {
 groups <- group_rows(d, keys)
 ans <- lapply(groups, function(ii) {
  if (length(ii) != 3L) stop("expected exactly three timing repeats for ", paste(d[ii[1],keys], collapse="/"), "; found ", length(ii))
  base <- d[ii[1],keys,drop=FALSE]
  do.call(rbind,lapply(names(metrics), function(metric) data.frame(base,metric=metric,mean_value=mean(d[[metrics[[metric]]]][ii]),count=length(ii),stringsAsFactors=FALSE)))
 })
 out <- do.call(rbind,ans); rownames(out) <- NULL; out
}
method_levels <- c("ordinary_rif_ols","gmoq","gmoq_jackknife")
method_labels <- c("ordinary RIF–OLS","G-MoQ pooled","G-MoQ jackknife")
method_cols <- c("ordinary RIF–OLS"="#1b9e77","G-MoQ pooled"="#d95f02","G-MoQ jackknife"="#7570b3","Mixture ML"="#1f78b4")
panel2 <- function(a, b) {
 panel_levels <- as.vector(t(outer(levels(a), levels(b), paste, sep=" | ")))
 factor(paste(a, b, sep=" | "), levels=panel_levels)
}
panel3 <- function(a, b, c) {
 panel_levels <- unlist(lapply(levels(a), function(aa) unlist(lapply(levels(b), function(bb) paste(aa, bb, levels(c), sep="\n")))))
 factor(paste(a, b, c, sep="\n"), levels=panel_levels)
}
missing <- character(); sources <- character()
# S1/S2: the reviewed N=400,000, M=1,000 Stage 1 benchmark.
s1d <- rootcsv(".local/stage1-figures-20260917-05/tables")
pairs <- num(readc(file.path(s1d,"S1_pairs.csv"),c("snp_id","tau","beta_gmoq","beta_rif_ols","difference")),c("tau","beta_gmoq","beta_rif_ols","difference"))
se <- num(readc(file.path(s1d,"S1_SE.csv"),c("snp_id","tau","se_gmoq","se_ordinary_rif_ols","se_gmoq_jackknife")),c("tau","se_gmoq","se_ordinary_rif_ols","se_gmoq_jackknife"))
tm <- num(readc(file.path(s1d,"S2_timings.csv"),c("workload","method","repeat","n","m","fit_seconds","wall_seconds","peak_rss_kib")),c("repeat","n","m","fit_seconds","wall_seconds","peak_rss_kib"))
if (length(unique(pairs$tau)) != 45L || nrow(tm) != 36L) stop("unexpected S1/S2 dimensions")
s2_cost <- unlist(tm[c("fit_seconds","wall_seconds","peak_rss_kib")], use.names=FALSE)
if (any(s2_cost <= 0)) stop("S2 timing/memory costs must be positive for log10 axes; zero-cost rows are not silently dropped")
copycsv(file.path(s1d,"S1_pairs.csv"),"S1_pairs.csv"); copycsv(file.path(s1d,"S1_SE.csv"),"S1_SE.csv"); copycsv(file.path(s1d,"S2_timings.csv"),"S2_timings.csv"); sources <- c(sources,s1d)
p <- ggplot(pairs,aes(beta_rif_ols,beta_gmoq,colour=tau))+geom_point(size=.55,alpha=.45)+geom_abline(slope=1,linetype=2)+scale_colour_gradient(low="#2166ac",high="#b2182b")+labs(title="S1. Stage 1 coefficient agreement",x="Ordinary RIF-OLS coefficient",y="G-MoQ coefficient",colour="tau")+plain
savep(p,"S1_beta",8,5.5)
sl <- rbind(data.frame(reference=se$se_ordinary_rif_ols,estimate=se$se_gmoq,source="Pooled"),data.frame(reference=se$se_ordinary_rif_ols,estimate=se$se_gmoq_jackknife,source="Jackknife"))
p <- ggplot(sl,aes(reference,estimate))+geom_point(size=.35,alpha=.25)+geom_abline(slope=1,linetype=2)+facet_wrap(~source,scales="free")+labs(title="S1. Stage 1 standard-error agreement",x="Ordinary RIF-OLS SE",y="G-MoQ SE")+plain
savep(p,"S1_SE",8,4.5)
tm$method <- factor(tm$method,levels=method_levels)
wi <- unique(tm[c("workload","n","m")]); wi <- wi[order(wi$n,wi$m),]; tm$workload <- factor(tm$workload,levels=wi$workload)
tm <- tm[order(tm$workload,tm$method,tm[["repeat"]]),,drop=FALSE]
s2means <- mean_timings(tm,c("workload","method","n","m"),c("fit + SE seconds"="fit_seconds","whole-process seconds"="wall_seconds","peak RSS MiB"="peak_rss_kib"))
s2means$mean_value[s2means$metric=="peak RSS MiB"] <- s2means$mean_value[s2means$metric=="peak RSS MiB"]/1024
if (any(s2means$count != 3L) || nrow(s2means) != 36L) stop("unexpected S2 mean timing dimensions")
write.csv(s2means,file.path(out,"S2_mean_timings.csv"),row.names=FALSE)
s2plot <- s2means; s2plot$method_label <- factor(as.character(s2plot$method),levels=method_levels,labels=method_labels); s2plot$workload <- factor(s2plot$workload,levels=levels(tm$workload))
p <- ggplot(s2plot,aes(workload,mean_value,colour=method_label))+geom_point(position=position_dodge(width=.7),size=2.2)+facet_wrap(~metric,scales="free",nrow=1)+scale_y_log10(labels=scales::label_number())+scale_colour_manual(values=method_cols,limits=method_labels,drop=FALSE)+labs(title="S2. Matched Stage 1 timing and memory",x="Workload",y="Mean cost (log10 scale; labels in original units)",colour="Method",caption="Each point is the arithmetic mean of three timing repeats; axes use log10 scaling with original-unit labels.")+plain+theme(axis.text.x=element_text(angle=35,hjust=1))
savep(p,"S2",10,4.8)
# S3: use the largest completed recovery bank; fall back only to a complete N=50,000 bank.
s3dirs <- c(rootcsv(".local/stage2-recovery-n400000-m1000-20260917-01"),rootcsv(".local/stage2-recovery-n50000-m1000-20260917-01")); s3d <- NULL
for (d in s3dirs) if (file.exists(file.path(d,"summary.json"))) { q <- fromJSON(file.path(d,"summary.json")); if (identical(q$status,"ok") && as.integer(q$attempted)==4000L && as.integer(q$failed)==0L) { s3d <- d; break } }
if (is.null(s3d)) missing <- c(missing,"S3: no complete 4,000-dataset recovery bank") else {
 s3 <- readc(file.path(s3d,"results.csv"),c("scenario","n","parameter","covariance_source","generating_effect","projection_target","estimate","se","ci_lower","ci_upper","status")); s3 <- num(s3,c("n","generating_effect","projection_target","estimate","se","ci_lower","ci_upper")); s3 <- s3[s3$status %in% c("success","ok"),]; active <- (s3$scenario=="normal_location"&s3$parameter=="beta_mu")|(s3$scenario=="normal_dispersion"&s3$parameter=="beta_sigma2")|(s3$scenario=="normal_both"&s3$parameter%in%c("beta_mu","beta_sigma2"))|(s3$scenario=="mixture_left"&s3$parameter=="beta_1"); s3 <- s3[active,]; copycsv(file.path(s3d,"results.csv"),"S3_results.csv"); sources <- c(sources,s3d)
 z <- rbind(data.frame(s3,truth=s3$generating_effect,truth_kind="Generating scientific effect"),data.frame(s3,truth=s3$projection_target,truth_kind="Fixed projection target")); z$scenario <- factor(z$scenario,levels=c("normal_location","normal_dispersion","normal_both","mixture_left"),labels=c("Mean effect","Variance effect","Mean + variance","Mixture: left mean")); z$covariance_source <- factor(z$covariance_source,levels=c("pooled","jackknife"),labels=c("Pooled","Jackknife")); z$panel <- panel2(z$covariance_source,z$scenario)
 p <- ggplot(z,aes(truth,estimate,colour=truth_kind,shape=parameter))+geom_point(size=.7,alpha=.45)+geom_abline(slope=1,linetype=2)+facet_wrap(~panel,scales="free",nrow=2)+labs(title="S3. Stage 2 recovery against two defined truths",subtitle=paste0("Complete bank; N=",format(unique(s3$n),big.mark=",",scientific=FALSE),"; failed fits are retained in source CSV"),x="Truth shown on x-axis",y="G-MoQ estimate",colour="Truth",shape="Parameter")+plain+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="bottom"); savep(p,"S3",12,7)
}
# S4/S5 are accepted only after the full mixture runner writes COMPLETE.json and final tables.
md <- rootcsv(".local/mixture-full-20260918-01"); mixok <- FALSE
if (file.exists(file.path(md,"COMPLETE.json"))) { cq <- fromJSON(file.path(md,"COMPLETE.json")); mixok <- identical(cq$status,"data_complete") && file.exists(file.path(md,"results.csv")) && file.exists(file.path(md,"panel_timings.csv")) }
if (!mixok) missing <- c(missing,"S4/S5: mixture-full final results/panel timings are not complete") else {
 mx <- readc(file.path(md,"results.csv"),c("dataset_id","n","replicate","method","parameter","estimate","se","known_truth","status")); mx <- num(mx,c("n","replicate","estimate","se","known_truth")); mx <- mx[mx$status=="ok",]; mx$n_label <- factor(mx$n,levels=sort(unique(mx$n)),labels=format(sort(unique(mx$n)),big.mark=",",scientific=FALSE)); copycsv(file.path(md,"results.csv"),"S4_results.csv"); copycsv(file.path(md,"panel_timings.csv"),"S5_panel_timings.csv"); sources <- c(sources,md); g <- mx[mx$method=="gmoq_pooled"&mx$parameter=="delta1",]; m <- mx[mx$method=="mixture_ml"&mx$parameter=="delta1",]; k <- merge(g,m,by=c("dataset_id","n","replicate","parameter"),suffixes=c("_gmoq","_ml")); k$n_label <- factor(k$n,levels=sort(unique(k$n)),labels=format(sort(unique(k$n)),big.mark=",",scientific=FALSE)); p <- ggplot(k,aes(estimate_ml,estimate_gmoq,colour=n_label))+geom_point(size=1,alpha=.7)+geom_abline(slope=1,linetype=2)+facet_wrap(~n_label,scales="free")+labs(title="S4. G-MoQ versus full mixture likelihood",x="Mixture ML estimate: left-component effect",y="G-MoQ pooled estimate",colour="N")+plain; savep(p,"S4_beta",9,4.8)
 t <- rbind(data.frame(estimate=g$estimate,truth=g$known_truth,method="G-MoQ pooled",n=g$n,n_label=g$n_label),data.frame(estimate=m$estimate,truth=m$known_truth,method="Mixture ML",n=m$n,n_label=m$n_label)); p <- ggplot(t,aes(truth,estimate,colour=method))+geom_point(size=1,alpha=.7)+geom_abline(slope=1,linetype=2)+facet_wrap(~n_label,scales="free")+labs(title="S4. Recovery against method-specific known truth",x="Known effect target",y="Estimate",colour="Method")+plain; savep(p,"S4_truth",9,4.8)
 p <- ggplot(k,aes(se_ml,se_gmoq,colour=n_label))+geom_point(size=1,alpha=.7)+geom_abline(slope=1,linetype=2)+facet_wrap(~n_label,scales="free")+labs(title="S4. Standard-error comparison",x="Mixture ML SE",y="G-MoQ pooled SE",colour="N")+plain; savep(p,"S4_SE",9,4.8)
 pt <- num(readc(file.path(md,"panel_timings.csv"),c("n","M","method","fit_seconds","wall_seconds","peak_rss_kib","attempted_snp_fits","successful_snp_fits","failures")),c("n","M","fit_seconds","wall_seconds","peak_rss_kib","attempted_snp_fits","successful_snp_fits","failures")); s5_cost <- unlist(pt[c("fit_seconds","wall_seconds","peak_rss_kib")], use.names=FALSE); if (any(s5_cost <= 0)) stop("S5 timing/memory costs must be positive for log10 axes; zero-cost rows are not silently dropped"); good <- pt$panel_status=="complete"&pt$successful_snp_fits==pt$attempted_snp_fits&pt$failures==0&pt$process_status=="success"; if (!any(good)) stop("mixture has no complete timing panels"); s5means <- mean_timings(pt[good,,drop=FALSE],c("n","M","method"),c("fit + SE seconds"="fit_seconds","whole-process seconds"="wall_seconds","peak RSS MiB"="peak_rss_kib")); s5means$mean_value[s5means$metric=="peak RSS MiB"] <- s5means$mean_value[s5means$metric=="peak RSS MiB"]/1024; if (any(s5means$count != 3L) || nrow(s5means) != 27L) stop("unexpected S5 mean timing dimensions"); write.csv(s5means,file.path(out,"S5_mean_timings.csv"),row.names=FALSE); tl <- s5means; tl$method <- factor(tl$method,levels=c("gmoq_pooled","gmoq_jackknife","mixture_ml"),labels=c("G-MoQ pooled","G-MoQ jackknife","Mixture ML")); tl$n_label <- factor(tl$n,levels=sort(unique(tl$n)),labels=format(sort(unique(tl$n)),big.mark=",",scientific=FALSE)); p <- ggplot(tl,aes(n_label,mean_value,colour=method))+geom_point(position=position_dodge(width=.7),size=2.2)+facet_wrap(~metric,scales="free",nrow=1)+scale_y_log10(labels=scales::label_number())+scale_colour_manual(values=method_cols,limits=c("G-MoQ pooled","G-MoQ jackknife","Mixture ML"),drop=FALSE)+labs(title="S5. Matched M=100 mixture-panel cost",subtitle="Only complete panels; independent-dataset timing is not the primary speed comparison",x="N",y="Mean cost (log10 scale; labels in original units)",colour="Method",caption="Each point is the arithmetic mean of three timing repeats; axes use log10 scaling with original-unit labels.")+plain; savep(p,"S5",10,4.8)
}
# S6: paired SE plots; the exact .02/.50/.98 Stage 1 taus are highlighted.
rd <- rootcsv(".local/rare-se-20260918-01"); r1 <- readc(file.path(rd,"stage1_tau.csv"),c("n","maf","scenario","tau","se_pooled","se_jackknife")); r1 <- num(r1,c("n","maf","tau","se_pooled","se_jackknife")); r1 <- r1[r1$scenario=="normal_null",]; r1$tau_group <- ifelse(r1$tau %in% c(.02,.98),"Extreme tau (.02/.98)",ifelse(r1$tau==.5,"Median tau (.50)","Other taus")); r1$n_label <- factor(r1$n,levels=c(10000,50000,400000),labels=c("N = 10,000","N = 50,000","N = 400,000")); r1$maf_label <- factor(r1$maf,levels=c(.05,.01,.001),labels=c("MAF = 5%","MAF = 1%","MAF = 0.1%")); r1$panel <- panel2(r1$n_label,r1$maf_label); copycsv(file.path(rd,"stage1_tau.csv"),"S6_stage1_tau.csv"); copycsv(file.path(rd,"stage1_summary.csv"),"S6_stage1_summary.csv"); sources <- c(sources,rd); p <- ggplot(r1,aes(se_pooled,se_jackknife,colour=tau_group))+geom_point(size=.45,alpha=.35)+geom_abline(slope=1,linetype=2)+facet_wrap(~panel,scales="free",nrow=3)+labs(title="S6A. Rare-variant Stage 1 SE sensitivity",subtitle="Null cells; extreme tails are explicitly shown",x="Pooled SE",y="Jackknife SE",colour=NULL)+plain; savep(p,"S6A",11,7)
r2 <- readc(file.path(rd,"stage2_parameters.csv"),c("n","maf","scenario","param","source","se")); r2 <- num(r2,c("n","maf","se")); r2 <- r2[r2$scenario%in%c("normal_null","normal_dispersion")&r2$source%in%c("pooled","jackknife"),]; key <- c("cell_index","replicate","n","maf","scenario","param"); po <- r2[r2$source=="pooled",c(key,"se")]; jk <- r2[r2$source=="jackknife",c(key,"se")]; r2p <- merge(po,jk,by=key,suffixes=c("_pooled","_jackknife")); r2p$n_label <- factor(r2p$n,levels=c(10000,50000,400000),labels=c("N = 10,000","N = 50,000","N = 400,000")); r2p$maf_label <- factor(r2p$maf,levels=c(.05,.01,.001),labels=c("MAF = 5%","MAF = 1%","MAF = 0.1%")); r2p$panel <- panel2(r2p$n_label,r2p$maf_label); p <- ggplot(r2p,aes(se_pooled,se_jackknife,colour=interaction(scenario,param,drop=TRUE)))+geom_point(size=.7,alpha=.45)+geom_abline(slope=1,linetype=2)+facet_wrap(~panel,scales="free",nrow=3)+labs(title="S6B. Rare-variant Stage 2 SE sensitivity",subtitle="Paired pooled and jackknife fits; coverage is reported in Table S2",x="Pooled SE",y="Jackknife SE",colour="Scenario / parameter")+plain; savep(p,"S6B",11,7); copycsv(file.path(rd,"stage2_parameters.csv"),"S6_stage2_parameters.csv"); copycsv(file.path(rd,"summary.csv"),"S6_summary.csv")
# S6C: estimate distributions use pooled rows only; OLS point estimates equal JK.
s6c_cols <- c("cell_index","replicate","n","maf","scenario","param","source","estimate","se","projection_target","generating_parameter","generating_target","projection_bias","actual_mac")
r3 <- readc(file.path(rd,"stage2_parameters.csv"),s6c_cols)
r3 <- num(r3,c("cell_index","replicate","n","maf","estimate","se","projection_target","generating_target","projection_bias","actual_mac"))
r3 <- r3[r3$scenario%in%c("normal_null","normal_dispersion")&r3$source=="pooled",]
if (nrow(r3)!=7200L) stop("unexpected S6C pooled effect dimensions: ",nrow(r3))
r3$maf_label <- factor(r3$maf,levels=c(.05,.01,.001),labels=c("MAF = 5%","MAF = 1%","MAF = 0.1%"))
r3$n_label <- factor(r3$n,levels=c(10000,50000,400000),labels=c("N = 10,000","N = 50,000","N = 400,000"))
r3$scenario_label <- factor(r3$scenario,levels=c("normal_null","normal_dispersion"),labels=c("Null","Dispersion signal (c = 0.1)"))
r3$parameter_label <- factor(r3$param,levels=c("beta_mu","beta_sigma2"),labels=c("Location (beta_mu)","Dispersion (beta_sigma2)"))
r3$panel <- panel3(r3$n_label,r3$scenario_label,r3$parameter_label)
truth <- unique(r3[c("n_label","scenario_label","parameter_label","maf_label","panel","projection_target")])
p <- ggplot(r3,aes(maf_label,estimate))
p <- p + geom_boxplot(width=.58,outlier.shape=NA,fill="#d9d9d9",colour="#555555")
p <- p + geom_jitter(position=position_jitter(width=.12,height=0,seed=13),alpha=.22,size=.45)
p <- p + geom_point(data=truth,aes(maf_label,projection_target),inherit.aes=FALSE,shape=95,size=8,colour="#b2182b")
p <- p + facet_wrap(~panel,scales="free",nrow=3)
p <- p + labs(title="S6C. Rare-variant Stage 2 estimate distributions",subtitle="Pooled OLS rows only; red horizontal marks are known population projection targets",x=NULL,y="Estimated effect",caption="200 replicates per cell; jackknife point estimates are identical under OLS and are not plotted")
p <- p + plain + theme(axis.text.x=element_text(angle=35,hjust=1),strip.text=element_text(size=8))
savep(p,"S6C",12,9)
write.csv(r3[s6c_cols],file.path(out,"S6C_effects.csv"),row.names=FALSE)
copycsv(file.path(rd,"summary.csv"),"S6C_summary.csv")
# Compact Table S1.
`%||%` <- function(x,y) if (length(x) && !is.null(x)) x else y
pin <- "9c47e7ba56370e708810fff8c2c9ad3b723cbccf"; s1tab <- data.frame(study=paste0("S",1:6),DGP=c("Matched Stage 1 RIF–OLS","Matched Stage 1 timing","Normal/mixture single-effect recovery","Mixture projection vs ML","Matched mixture panel","Rare-variant null/dispersion"),N=c("400000","400000","400000 (fallback 50000)","10000/50000/400000","10000/50000/400000","10000/50000/400000"),M_MAF=c("M=1,1000","M=1/256/1000","MAF .30; M=1000","MAF .30; M=100 panel","M=100 panel","MAF .05/.01/.001"),replicates=c("3 timing repeats","3 timing repeats","1000/scenario","20/N","3 panels/N","200/cell"),effect_range=c("NA","NA","U(-.1,.1)","left component U(-.1,.1)","shared panel genotypes","null; c=.1"),tau45=c("yes","yes","yes","yes","no","yes Stage 1"),blocks1000=c("yes","yes","yes","yes","yes","yes"),package_pin=pin,target=c("coefficient equality","fit/whole-process/RSS","generating and projection","method-specific known truth","cost only","projection and paired SE"),threads=c("1","1","1","1","1","1"),source=c(s1d, s1d, s3d %||% "missing",md,md,rd),check.names=FALSE,stringsAsFactors=FALSE)
write.csv(s1tab,file.path(out,"TableS1_design_provenance.csv"),row.names=FALSE)
# Table S2 uses one tidy schema across numerical recovery, SE agreement, and timings.
one <- function(study,quantity,scenario=NA,n=NA,maf=NA,parameter=NA,method=NA,cov=NA,workload=NA,bias=NA,rmse=NA,esd=NA,mse=NA,covg=NA,mcse=NA,attempted=NA,success=NA,fail=NA,seagr=NA,fit=NA,wall=NA,rss=NA,timing_repeats=NA,source=NA) data.frame(study,quantity,scenario,n,maf,parameter,method,covariance_source=cov,workload,bias,rmse,empirical_sd=esd,mean_se=mse,coverage=covg,coverage_mcse=mcse,attempted,success,fail,se_agreement=seagr,timing_fit_seconds=fit,timing_wall_seconds=wall,peak_rss_mib=rss,timing_repeats,source,stringsAsFactors=FALSE)
tab <- list();
for (ii in group_rows(tm,c("workload","method"))) { q <- tm[ii,,drop=FALSE]; if (nrow(q) != 3L) stop("unexpected S2 timing repeats in Table S2"); tab[[length(tab)+1L]] <- one("S2","matched_stage1_timing_mean",n=q$n[1],method=as.character(q$method[1]),workload=as.character(q$workload[1]),attempted=3L,success=3L,fail=0L,fit=mean(q$fit_seconds),wall=mean(q$wall_seconds),rss=mean(q$peak_rss_kib)/1024,timing_repeats=3L,source=s1d) }
tab[[length(tab)+1L]] <- one("S1","coefficient_minus_RIF_OLS",bias=mean(pairs$difference),rmse=sqrt(mean(pairs$difference^2)),esd=sd(pairs$difference),attempted=nrow(pairs),success=nrow(pairs),source=s1d); for (v in c("se_gmoq","se_gmoq_jackknife")) { rr <- se[[v]]/se$se_ordinary_rif_ols; tab[[length(tab)+1]] <- one("S1","SE_ratio_to_ordinary",method=ifelse(v=="se_gmoq","pooled","jackknife"),bias=mean(rr-1),rmse=sqrt(mean((rr-1)^2)),esd=sd(rr),seagr=median(rr),attempted=length(rr),success=length(rr),source=s1d) }
if (!is.null(s3d)) { for (i in split(seq_len(nrow(s3)),interaction(s3$scenario,s3$covariance_source,s3$parameter,drop=TRUE))) { q <- s3[i,]; ok <- is.finite(q$estimate); e <- q$estimate[ok]-q$projection_target[ok]; cv <- mean(q$ci_lower[ok]<=q$projection_target[ok]&q$ci_upper[ok]>=q$projection_target[ok]); tab[[length(tab)+1]] <- one("S3","estimate_minus_projection",q$scenario[1],q$n[1],parameter=q$parameter[1],method="G-MoQ",cov=q$covariance_source[1],bias=mean(e),rmse=sqrt(mean(e^2)),esd=sd(q$estimate[ok]),mse=mean(q$se[ok]),covg=cv,mcse=sqrt(cv*(1-cv)/sum(ok)),attempted=nrow(q),success=sum(ok),fail=sum(!ok),source=s3d) } }
if (mixok) { for (i in split(seq_len(nrow(mx)),interaction(mx$n,mx$method,mx$parameter,drop=TRUE))) { q <- mx[i,]; e <- q$estimate-q$known_truth; tab[[length(tab)+1]] <- one("S4","estimate_minus_method_truth",method=q$method[1],n=q$n[1],parameter=q$parameter[1],bias=mean(e),rmse=sqrt(mean(e^2)),esd=sd(q$estimate),mse=mean(q$se),attempted=nrow(q),success=nrow(q),source=md) } }
if (mixok) { for (ii in group_rows(pt,c("n","M","method"))) { q <- pt[ii,,drop=FALSE]; tab[[length(tab)+1]] <- one("S5","matched_panel_timing_mean",method=q$method[1],n=q$n[1],workload=paste0("M=",q$M[1]),fit=mean(q$fit_seconds,na.rm=TRUE),wall=mean(q$wall_seconds,na.rm=TRUE),rss=mean(q$peak_rss_kib,na.rm=TRUE)/1024,attempted=sum(q$attempted_snp_fits),success=sum(q$successful_snp_fits),fail=sum(q$failures),timing_repeats=nrow(q),source=md) } }
s6sum <- readc(file.path(rd,"summary.csv")); for (i in seq_len(nrow(s6sum))) { q <- s6sum[i,]; tab[[length(tab)+1]] <- one("S6","estimate_minus_projection",scenario=q$scenario,n=as.numeric(q$n),maf=as.numeric(q$maf),parameter=q$param,method=q$source,bias=as.numeric(q$projectionbias),esd=as.numeric(q$empiricalSD),mse=as.numeric(q$meanSE),covg=as.numeric(q$coverage),mcse=as.numeric(q$coverageMCSE),attempted=as.numeric(q$attempted),success=as.numeric(q$success),fail=as.numeric(q$fail),seagr=as.numeric(q$medianpairedSEratio),source=rd) }
write.csv(do.call(rbind,tab),file.path(out,"TableS2_numerical_summary.csv"),row.names=FALSE)
# Provenance, captions, and a short explicit completion record.
commit <- tryCatch(system2("git",c("-C",root,"rev-parse","HEAD"),stdout=TRUE),error=function(e) "unavailable")
cap <- c("S1: Coefficient and standard-error agreement with independently fitted ordinary RIF–OLS; SE identity is descriptive because covariance conventions differ.","S2: Matched Stage 1 fit/SE, whole-process timing, and peak RSS; each point is the arithmetic mean of three timing repeats in original units, with log10 axes and original-unit labels.","S3: Estimated effects against generating scientific effects and fixed population projection targets; their gap is not estimation bias.","S4: G-MoQ pooled and full mixture likelihood on the same datasets, with method-specific known truths; all parameters remain in Table S2.","S5: Complete matched M=100 panel timing only; each point is the arithmetic mean of three timing repeats in original units, with fit+SE, whole-process wall time, and RSS shown on log10 axes with original-unit labels.","S6A/B: Paired pooled versus jackknife SEs across N and MAF. Extreme Stage 1 taus (.02/.98) and median tau (.50) are highlighted; neither covariance source is automatically calibrated.","S6C: Distributions of pooled OLS Stage 2 estimates across 200 replicates per cell for the normal null and c = 0.1 dispersion signal. Red horizontal marks show the known population projection target for each cell, which varies slightly by MAF; all N, MAF and both parameters are shown, and jackknife rows are excluded because OLS point estimates are identical.","Faceted panels use independently scaled axes; coefficient and truth-target panels remain linear, while S2/S5 timing and memory cost axes use log10 scaling with original-unit tick labels.")
if (length(missing)) cap <- c(cap,"Not rendered because inputs were incomplete:",paste0("- ",missing))
writeLines(cap,file.path(out,"captions.txt")); writeLines(c("Simulation supplement render",paste0("output=",out),paste0("paper_commit=",commit),paste0("package_pin=",pin),paste0("sources=",paste(unique(sources),collapse=";")),paste0("underlying_csvs=",canonical," (cluster copy; large source CSVs remain on the cluster)"),if(length(missing)) c("missing=",missing) else "completion=all approved figures rendered","All plotted values are from aggregate CSVs; no participant-level data are exported."),file.path(out,"README.md")); writeLines(c(paste0("R=",R.version.string),paste0("ggplot2=",packageVersion("ggplot2")),paste0("root=",root),paste0("output=",out),paste0("underlying_csvs=",canonical)),file.path(out,"render_metadata.txt"))
