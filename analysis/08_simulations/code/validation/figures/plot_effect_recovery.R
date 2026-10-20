#!/usr/bin/env Rscript
# Draft plots of Stage 2 effect recovery
#
# What: plots estimates against the generating effect and against the projection target for
#   the study 3 runs, as a run-time check. The paper's Supplementary Fig. S3 is drawn by
#   figures/figS03_projection_targets.R from the imported summaries.
# Paper: none directly (study 3)
# Input: results.csv written by simulate_stage2.py
# Output: working plots and effect_recovery_summary.csv in a new directory
# Runs: on a cluster, called by stage2.sbatch, under a minute
# Usage: Rscript plot_effect_recovery.R results.csv NEW_OUTPUT_DIR
args <- commandArgs(TRUE)
stopifnot(length(args)==2L, !file.exists(args[2]))
suppressPackageStartupMessages(library(ggplot2))
dir.create(args[2], recursive=TRUE)
d <- read.csv(args[1], stringsAsFactors=FALSE)
needed <- c("scenario","snp_id","parameter","covariance_source","generating_effect",
            "projection_target","estimate","se","ci_lower","ci_upper","status")
stopifnot(all(needed %in% names(d)))
groups <- split(d, interaction(d$scenario,d$covariance_source,d$parameter,drop=TRUE))
summary <- do.call(rbind,lapply(groups,function(g) {
  z <- g[g$status %in% c("ok","success"),]
  e <- z$estimate-z$projection_target
  data.frame(scenario=g$scenario[1],covariance_source=g$covariance_source[1],parameter=g$parameter[1],
    n_attempted=nrow(g),n_success=nrow(z),n_failed=nrow(g)-nrow(z),
    projection_bias=mean(e),projection_rmse=sqrt(mean(e^2)),
    bias_mcse=sd(e)/sqrt(nrow(z)),generating_bias=mean(z$estimate-z$generating_effect),
    generating_rmse=sqrt(mean((z$estimate-z$generating_effect)^2)),mean_se=mean(z$se),
    projection_coverage=mean(z$ci_lower<=z$projection_target & z$ci_upper>=z$projection_target))
}))
write.csv(summary,file.path(args[2],"effect_recovery_summary.csv"),row.names=FALSE)
z <- d[d$status %in% c("ok","success"),]
active <- with(z, (scenario=="normal_location" & parameter=="beta_mu") |
  (scenario=="normal_dispersion" & parameter=="beta_sigma2") |
  (scenario=="normal_both" & parameter %in% c("beta_mu","beta_sigma2")) |
  (scenario=="mixture_left" & parameter=="beta_1"))
z <- z[active,]
z$scenario <- factor(z$scenario,levels=c("normal_location","normal_dispersion","normal_both","mixture_left"),
                     labels=c("Mean effect","Variance effect","Mean + variance","Mixture: lower-component mean"))
for (source in c("pooled","jackknife")) {
  a <- z[z$covariance_source==source,]
  stopifnot(nrow(a)>0,all(is.finite(a$estimate)),all(is.finite(a$projection_target)))
  for (target in c("generating_effect","projection_target")) {
    a$truth <- a[[target]]
    p <- ggplot(a,aes(truth,estimate,colour=parameter)) +
      geom_abline(slope=1,intercept=0,linetype=2,colour="grey40") +
      geom_point(size=1.6,alpha=.65) + facet_wrap(~scenario,ncol=2) +
      coord_equal() + theme_bw(base_size=11) + theme(legend.position="bottom") +
      labs(title=paste("Stage 2 effect recovery:",source,"OLS"),
        subtitle=paste0("N=",format(unique(a$n),big.mark=","),"; independent datasets across SNPs"),
        x=if(target=="generating_effect") "Generating scientific effect" else "Population projection target",
        y="G-MoQ estimate",colour="Parameter")
    stem <- paste0("S3_",source,"_",target)
    ggsave(file.path(args[2],paste0(stem,".pdf")),p,width=8,height=7)
    ggsave(file.path(args[2],paste0(stem,".png")),p,width=8,height=7,dpi=300)
  }
}
writeLines(c("Generating effects and population projection targets are shown separately.",
 "Mean/variance/mixture axes are not assumed to equal generating effects exactly.",
 "Inactive parameters and failed fits are retained in the source CSV and summary.",
 paste0("R=",R.version.string,"; ggplot2=",packageVersion("ggplot2"))),file.path(args[2],"caption.txt"))
