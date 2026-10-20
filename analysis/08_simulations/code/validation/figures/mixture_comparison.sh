#!/usr/bin/env bash
# Mixture recovery and cost: G-MoQ against direct mixture maximum likelihood
#
# What: two-component normal mixture (equal weights, component means -2 and 2, unit SDs),
#   genotype MAF 0.3, genotype effect on the lower component location (delta1) drawn from
#   U(-0.1, 0.1). For each of 20 data sets at N = 10,000, 50,000 and 400,000, fits the
#   membership and both component locations three ways, each in a fresh process timed by
#   GNU time: G-MoQ (package Stage 1, then Stage 2 with the analytic mixture weights) with the
#   pooled covariance, the same with the 1,000-block jackknife, and the conditional mixture
#   likelihood refitted by maximum likelihood (comparators.fit_mixture_ml, observed-Hessian
#   SEs). G-MoQ is compared with its exact population projection target (truth.py), maximum
#   likelihood with the generating effect. Then, for the cost comparison, three matched
#   panels per N of 100 null variants sharing one phenotype are fitted whole by each method:
#   G-MoQ runs Stage 1 once for the panel and Stage 2 in one batch (param_gwas), maximum
#   likelihood refits the mixture variant by variant.
# Paper: Supplementary Figs. S4a,b,d and S5a-c; Supplementary Table S2 (studies 4 and 5)
# Input: none (simulated data); ../comparators.py, ../truth.py, ../run_examples.py,
#   ../stage2.R; the G-MoQ packages through the environment made by ../bootstrap.sh
# Output: <ROOT>/.local/<run id>/{results.csv,summary.csv,failures.csv,panel_timings.csv,
#   panel_estimates.csv,provenance.json,COMPLETE.json,figures/}
# Runs: on a cluster, one core, about 2.5 hours:
#   sbatch validation/figures/mixture_comparison.sh mixture-full-20260918-01 20 \
#     10000,50000,400000 --panel-size 100 --panel-reps 3
#   Seeds: data set i (from 0) uses seed 96,000,000 + i and effect seed 95,000,000 + i (the
#   same data sets at every N); panel i uses phenotype seed 97,000,000 + i and genotype
#   seed 97,000,000 + i + 7,919.
#   The same file is re-entered as a child process for each timed arm (--child,
#   --panel-child), for the batch Stage 2 of a panel (--panel-stage2, embedded R) and for
#   working plots (--plot). The #SBATCH --chdir and --output lines hold the authors' cluster
#   paths; edit them, and the user settings below, for your cluster.
#SBATCH --job-name=gmoq-mixture
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --chdir=/user/work/fh6520/Gmoq_paper-figures-20260917
#SBATCH --output=/user/work/fh6520/Gmoq_paper-figures-20260917/.local/mixture-%j.out
set -Eeuo pipefail
## User settings
# Directory that contains validation/ (here: analysis/08_simulations/code).
ROOT=${GMOQ_SIM_ROOT:-/user/work/fh6520/Gmoq_paper-figures-20260917}
# Environment made by validation/bootstrap.sh.
BOOT=${GMOQ_SIM_ENV:-/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-bootstrap-20260917-02}
export GMOQ_SIM_ROOT="$ROOT" GMOQ_SIM_ENV="$BOOT"
source "$BOOT/activate.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VALIDATION_COMPUTE_HOST=1 MIXTURE_SCRIPT_PATH="$ROOT/validation/figures/mixture_comparison.sh"
if [[ "${1:-}" == --panel-stage2 ]]; then
  shift
  exec "$GMOQ_RSCRIPT_BIN" --vanilla - "$@" <<'RS2'
a <- commandArgs(trailingOnly=TRUE); k <- a[seq(1,length(a),2)]; v <- a[seq(2,length(a),2)]
z <- as.list(v); names(z) <- sub("^--","",k)
need <- function(n) { if (is.null(z[[n]]) || !nzchar(z[[n]])) stop("missing --",n); z[[n]] }
if (!requireNamespace("jsonlite",quietly=TRUE)) stop("jsonlite is required")
root <- normalizePath(need("package-root"),mustWork=TRUE); got <- trimws(system2("git",c("-C",root,"rev-parse","HEAD"),stdout=TRUE))
if (!identical(got,need("package-commit"))) stop("panel Stage 2 source pin mismatch")
rlib <- normalizePath(need("r-lib"),mustWork=TRUE); .libPaths(unique(c(rlib,.libPaths())))
library(gmoq,lib.loc=rlib); if (as.character(packageVersion("gmoq"))!="0.4.0") stop("wrong gmoq version")
x <- jsonlite::fromJSON(need("input"),simplifyVector=TRUE); source <- need("source")
stage <- list(taus=as.numeric(x$taus),q_tau=as.numeric(x$q_tau),Q_slope=as.matrix(x$beta),SE_tau=as.matrix(x$se_tau))
b <- x$baseline; W <- gmoq::make_weights_normal_mixture(stage$taus,stage$q_tau,as.numeric(b$p1),as.numeric(b$mu1),as.numeric(b$sd1),as.numeric(b$mu2),as.numeric(b$sd2),include_membership=TRUE)
cv <- as.numeric(x$cov_vec); off <- if (source=="pooled") rep(0,ncol(stage$Q_slope)) else as.numeric(x$offsets); sc <- if (source=="pooled") as.numeric(x$cov_scales) else NULL
nb <- if (source=="jackknife") 1000L else NULL; nu <- if (source=="pooled") as.numeric(x$n) else NULL
fit <- gmoq::param_gwas(stage1=stage,transform="custom_W",transform_args=list(W=W),se_mode="tau_cov",cov_data=list(cov_vec=cv,offsets=off,cov_scales=sc),tau_cov_estimator="ols",n_threads=1L,return_cov="none",q_method="saddlepoint",nu=nu,nu_source=if(source=="pooled") "panel_pooled" else NULL,n_blocks=nb,n_eff=as.numeric(x$n_eff),unknown_nu="error")
if (any(!is.finite(fit$params)) || any(!is.finite(fit$SE_params)) || any(fit$SE_params<=0)) stop("non-finite panel Stage 2 estimate/SE")
rows <- lapply(seq_len(ncol(fit$params)),function(j) list(snp_index=j-1L,estimates=as.numeric(fit$params[,j]),ses=as.numeric(fit$SE_params[,j])))
jsonlite::write_json(list(status="ok",parameter_names=rownames(fit$params),rows=rows),need("output"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
RS2
fi
if [[ "${1:-}" == --plot ]]; then
  shift
  exec "$GMOQ_RSCRIPT_BIN" --vanilla - "$@" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2L || length(args)>3L) stop("plot needs results.csv, figure directory, and optional panel_timings.csv")
if (!requireNamespace("ggplot2",quietly=TRUE)) stop("ggplot2 is required")
suppressPackageStartupMessages(library(ggplot2)); x <- read.csv(args[1],stringsAsFactors=FALSE,check.names=FALSE)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
x$method <- factor(x$method,levels=c("gmoq_pooled","gmoq_jackknife","mixture_ml"),labels=c("G-MoQ pooled","G-MoQ jackknife","conditional mixture ML"))
ok <- x[x$status=="ok" & is.finite(x$estimate) & is.finite(x$known_truth) & is.finite(x$se),]
savep <- function(p,n,w=9,h=6) { ggsave(file.path(args[2],paste0(n,".pdf")),p,width=w,height=h); ggsave(file.path(args[2],paste0(n,".png")),p,width=w,height=h,dpi=300) }
p <- ggplot(ok,aes(known_truth,estimate,colour=method)) + geom_point(alpha=.65,size=1.1) + geom_abline(slope=1,intercept=0,linetype=2) + facet_grid(parameter~n,scales="free") + theme_bw() + labs(title="Mixture recovery against the method-specific known target",subtitle="G-MoQ is assessed against its finite-grid population projection; ML against the generating scientific effect",x="Known target",y="Estimate",colour="Method")
savep(p,"S4_scatter",10,7)
u <- unique(ok[c("dataset_id","n","method","fit_seconds","wall_seconds","peak_rss_kib")])
l <- rbind(data.frame(u,metric="fit + SE seconds",value=u$fit_seconds),data.frame(u,metric="full-process wall seconds",value=u$wall_seconds),data.frame(u,metric="peak RSS MiB",value=u$peak_rss_kib/1024))
p <- ggplot(l,aes(method,value,colour=method)) + geom_point(position=position_jitter(width=.08),alpha=.7) + facet_grid(metric~n,scales="free_y") + theme_bw() + theme(axis.text.x=element_text(angle=25,hjust=1)) + labs(title="Independent-dataset mixture comparison cost",subtitle="Sensitivity: one fresh phenotype and refit per SNP; GNU time -v",x=NULL,y=NULL,colour="Method")
savep(p,"S5_timing_independent",10,7)
if (length(args)==3L && file.exists(args[3])) {
 pt <- read.csv(args[3],stringsAsFactors=FALSE,check.names=FALSE)
 if (nrow(pt)>0L) {
  pt$method <- factor(pt$method,levels=c("gmoq_pooled","gmoq_jackknife","mixture_ml"),labels=c("G-MoQ pooled","G-MoQ jackknife","conditional mixture ML"))
  pt <- pt[pt$panel_status=="complete" & pt$successful_snp_fits==pt$attempted_snp_fits & is.finite(pt$fit_seconds) & is.finite(pt$wall_seconds) & is.finite(pt$peak_rss_kib),]
  if (nrow(pt)>0L) {
   pl <- rbind(data.frame(pt,metric="fit + SE seconds",value=pt$fit_seconds),data.frame(pt,metric="full-process wall seconds",value=pt$wall_seconds),data.frame(pt,metric="peak RSS MiB",value=pt$peak_rss_kib/1024))
   pp <- ggplot(pl,aes(method,value,colour=method)) + geom_point(position=position_jitter(width=.08),size=1.8,alpha=.8) + facet_grid(metric~n,scales="free_y") + theme_bw() + theme(axis.text.x=element_text(angle=25,hjust=1)) + labs(title="Matched GWAS-panel mixture comparison cost",subtitle="Primary S5: one shared Y and M independent genotypes; failed panels are retained but excluded",x=NULL,y=NULL,colour="Method")
   savep(pp,"S5_timing",10,7)
  }
 }
}
writeLines(c("ggplot2 mixture comparison figures",paste0("R=",R.version.string),paste0("ggplot2=",packageVersion("ggplot2")),"S5 primary uses matched-panel timings; independent-dataset timing is S5_timing_independent and is sensitivity only.","Panel null results are timing/finite-SE checks, not calibration evidence."),file.path(args[2],"plot_metadata.txt"))
RS
fi
exec "$GMOQ_STAGE1_PYTHON" - "$@" <<'PY'
from __future__ import annotations
import argparse,csv,hashlib,json,os,pickle,signal,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
ROOT=Path(os.environ['GMOQ_SIM_ROOT']); sys.path.insert(0,str(ROOT))
from validation import truth
from validation.comparators import fit_mixture_ml
from validation.run_examples import (PINNED_PACKAGE_COMMIT,PINNED_PACKAGE_ROOT,_run_stage2,
 balanced_block_ids,load_environment_manifest,load_pinned_core)
TAUS=np.linspace(.02,.98,45); ENV=Path(os.environ['GMOQ_SIM_ENV'])/'environment.json'; METHODS=("gmoq_pooled","gmoq_jackknife","mixture_ml"); PARAMS=("gamma","delta1","delta2")

def safe(x):
 if isinstance(x,dict): return {str(k):safe(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)): return [safe(v) for v in x]
 if isinstance(x,np.ndarray): return safe(x.tolist())
 if isinstance(x,(np.integer,)): return int(x)
 if isinstance(x,(np.floating,)): return float(x)
 return x

def write(p,x):
 p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(safe(x),indent=2,allow_nan=False)+"\n")
def sha(p):
 h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
def git(path,args=("rev-parse","HEAD")): return subprocess.check_output(["git","-C",str(path),*args],text=True).strip()
def canon(name): return {"beta_1":"delta1","beta_2":"delta2","delta_1":"delta1","delta_2":"delta2"}.get(name,name)
def elapsed(s):
 z=s.strip().split(":")
 try: return sum(float(v)*60**(len(z)-i-1) for i,v in enumerate(z))
 except ValueError: return None
def timeinfo(p):
 d={}
 if not p.is_file(): return d
 for line in p.read_text(errors="replace").splitlines():
  if "Maximum resident set size" in line:
   try: d["peak_rss_kib"]=int(line.rsplit(":",1)[1].strip())
   except ValueError: pass
  if "Elapsed (wall clock) time" in line:
   try: d["wall_seconds"]=elapsed(line.split(": ",1)[1])
   except (IndexError,ValueError): pass
 return d

def child(a):
 payload=pickle.loads(Path(a.input).read_bytes()); g=np.asarray(payload["g"],float); y=np.asarray(payload["y"],float); pop=payload["truth"]; started=time.perf_counter(); arm=Path(a.workdir); arm.mkdir(parents=True,exist_ok=True)
 try:
  if a.method=="mixture_ml":
   fit=fit_mixture_ml(g,y,{"family":"mixture","maf":.3}); dt=time.perf_counter()-started; status=str(fit["status"]); est=fit.get("estimates",{}); ses=fit.get("standard_errors",{}); rows=[]
   for p in PARAMS: rows.append({"parameter":p,"estimate":float(est[p]) if status=="ok" and p in est else None,"se":float(ses[p]) if status=="ok" and p in ses else None,"status":status,"reason":fit.get("reason")})
   out={"status":"ok" if status=="ok" else status,"rows":rows,"fit_seconds":dt,"cost_label":"conditional mixture ML fit plus observed-Hessian SE"}
  else:
   env=load_environment_manifest(ENV,Path(PINNED_PACKAGE_ROOT)); core,prov,sigma0fn=load_pinned_core(Path(PINNED_PACKAGE_ROOT),env)
   base=truth.fit_baseline(y,"mixture"); rif,q,fhat=core.compute_rif(y,TAUS,density_method="silverman_gaussian",density_floor=1e-8)
   Qrif,Q=core.residualize(rif,np.ones((len(y),1))); gg=g[:,None]
   beta,d=core.compute_pooled_panel(gg,Qrif,Q,min_mac=5.,min_var=1e-8,n_threads=1,panel_size=16); beta=np.asarray(beta,float)[0]; d=float(np.asarray(d)[0])
   if a.method=="gmoq_pooled":
    sigma=np.asarray(sigma0fn(Qrif),float); se=np.sqrt(np.maximum(np.diag(sigma)/d,0.)); source="pooled"; cov={"sigma0":sigma}
   else:
    blocks=1000; bid=balanced_block_ids(len(y),blocks,int(a.seed)); plan=core.make_block_plan(Qrif,Q,bid,blocks,n_threads=1); jb,js,ju=core.compute_block_finalize_panel_fast(gg,Qrif,Q,bid,blocks,min_mac=5.,min_var=1e-8,n_threads=1,panel_width=16,plan=plan)
    beta=np.asarray(jb,float)[0]; se=np.asarray(js,float)[0]; source="jackknife"; cov={"cov_upper_jackknife":np.asarray(ju,float)[0]}
   payload2={"case":{"family":"mixture","maf":.3},"n":len(y),"seed":int(a.seed),"family":"mixture","taus":TAUS,"q_tau":q,"fhat_tau":fhat,"baseline_fit":base,"population_truth":pop,"beta_pooled":beta if source=="pooled" else np.zeros(45),"se_tau_pooled":se if source=="pooled" else np.ones(45),"beta_jackknife":beta if source=="jackknife" else np.zeros(45),"se_tau_jackknife":se if source=="jackknife" else np.ones(45),"d_resid":d,"n_blocks":1000,"density_method":"silverman_gaussian","density_convention":"pinned public Stage 1 compute_rif","stage2_weight_convention":"pinned Stage 2 analytic mixture weights","package":prov,**cov}
   inp=arm/"stage1_summary.json"; write(inp,payload2); outpath=arm/("stage2-"+source+".json"); fit=_run_stage2(env["rscript"],ROOT/"validation/stage2.R",inp,outpath,arm,source,Path(PINNED_PACKAGE_ROOT),PINNED_PACKAGE_COMMIT,"0.4.0",Path(env["path"]),env["r_lib"],1100.)
   vals=fit["estimate"]["estimate"]; errs=fit["estimate"]["se"]; names=fit["estimate"]["parameter_names"]; got={canon(n):(float(vals[n]),float(errs[n])) for n in names}; rows=[{"parameter":p,"estimate":got[p][0],"se":got[p][1],"status":"ok","reason":None} for p in PARAMS]
   out={"status":"ok","rows":rows,"fit_seconds":time.perf_counter()-started,"cost_label":"G-MoQ minimal public Stage 1 kernel plus pinned Stage 2 and child R"}
  print(json.dumps(safe(out),allow_nan=False)); return 0
 except Exception as e:
  dt=time.perf_counter()-started; print(json.dumps({"status":"fit_failure","reason":f"{type(e).__name__}: {e}","fit_seconds":dt,"rows":[{"parameter":p,"estimate":None,"se":None,"status":"fit_failure","reason":str(e)} for p in PARAMS]})); traceback.print_exc(); return 1

def arm(method,g,y,pop,seed,root):
 d=root/"runs"/f"n{len(y)}_rep{seed}"/method; d.mkdir(parents=True,exist_ok=False); tf=d/"time-v.txt"; so=d/"stdout.log"; se=d/"stderr.log"; payload=d.parent/"input.pkl"
 if not payload.exists(): payload.write_bytes(pickle.dumps({"g":g,"y":y,"truth":pop},protocol=pickle.HIGHEST_PROTOCOL))
 cmd=["/usr/bin/time","-v","-o",str(tf),"bash",os.environ["MIXTURE_SCRIPT_PATH"],"--child",method,"--workdir",str(d),"--seed",str(seed),"--input",str(payload)]
 try:
  with so.open("w") as fo,se.open("w") as fe:
   q=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=fo,stderr=fe,start_new_session=True)
   try: q.communicate(timeout=1200.); rc=q.returncode; status="success" if rc==0 else "failed"
   except subprocess.TimeoutExpired:
    os.killpg(q.pid,signal.SIGKILL); q.wait(); rc=None; status="timeout"
 except Exception:
  rc=None; status="launch_failure"
 raw={}
 try:
  lines=so.read_text().splitlines(); raw=json.loads(lines[-1]) if lines else {}
 except Exception: pass
 ti=timeinfo(tf); return raw,dict(ti,process_status=status,return_code=rc,arm_status=status,log_dir=str(d))

def panel_child(a):
 payload=pickle.loads(Path(a.input).read_bytes()); G=np.asarray(payload["G"],float); y=np.asarray(payload["y"],float); M=G.shape[1]; started=time.perf_counter(); arm=Path(a.workdir); arm.mkdir(parents=True,exist_ok=True)
 def bad(reason,status="failed"):
  return {"status":status,"attempted_snp_fits":M,"successful_snp_fits":0,"failures":M,"fit_seconds":time.perf_counter()-started,"rows":[{"snp_index":j,"parameter":p,"estimate":None,"se":None,"status":status,"reason":reason} for j in range(M) for p in PARAMS],"reason":reason}
 try:
  if a.method=="mixture_ml":
   rows=[]; errors=[]
   for j in range(M):
    fit=fit_mixture_ml(G[:,j],y,{"family":"mixture","maf":.3}); good=str(fit["status"])=="ok" and all(p in fit.get("estimates",{}) and p in fit.get("standard_errors",{}) and np.isfinite(fit["estimates"][p]) and np.isfinite(fit["standard_errors"][p]) and fit["standard_errors"][p]>0 for p in PARAMS)
    if not good: errors.append({"snp_index":j,"status":fit.get("status"),"reason":fit.get("reason")})
    for p in PARAMS: rows.append({"snp_index":j,"parameter":p,"estimate":float(fit.get("estimates",{}).get(p)) if good else None,"se":float(fit.get("standard_errors",{}).get(p)) if good else None,"status":"ok" if good else str(fit.get("status","fit_failure")),"reason":fit.get("reason")})
   write(arm/"panel_ml_errors.json",errors); good=M-len(errors); out={"status":"complete" if good==M else "failed","attempted_snp_fits":M,"successful_snp_fits":good,"failures":M-good,"fit_seconds":time.perf_counter()-started,"rows":rows,"reason":None if good==M else "one or more mixture-ML SNP fits failed"}
  else:
   env=load_environment_manifest(ENV,Path(PINNED_PACKAGE_ROOT)); core,prov,sigma0fn=load_pinned_core(Path(PINNED_PACKAGE_ROOT),env); base=truth.fit_baseline(y,"mixture"); rif,q,_=core.compute_rif(y,TAUS,density_method="silverman_gaussian",density_floor=1e-8); Qrif,Q=core.residualize(rif,np.ones((len(y),1))); beta,d=core.compute_pooled_panel(G,Qrif,Q,min_mac=5.,min_var=1e-8,n_threads=1,panel_size=16); beta=np.asarray(beta,float); d=np.asarray(d,float)
   if a.method=="gmoq_pooled":
    sigma=np.asarray(sigma0fn(Qrif),float); se=np.sqrt(np.maximum(np.diag(sigma)[None,:]/d[:,None],0.)); source="pooled"; cov={"cov_vec":sigma[np.triu_indices(45)],"offsets":np.zeros(M,dtype=int),"cov_scales":1./d}
   else:
    bid=balanced_block_ids(len(y),1000,int(a.seed)); plan=core.make_block_plan(Qrif,Q,bid,1000,n_threads=1); jb,js,ju=core.compute_block_finalize_panel_fast(G,Qrif,Q,bid,1000,min_mac=5.,min_var=1e-8,n_threads=1,panel_width=16,plan=plan); beta=np.asarray(jb,float); se=np.asarray(js,float); source="jackknife"; L=45*46//2; cov={"cov_vec":np.asarray(ju,float).reshape(-1) ,"offsets":np.arange(M,dtype=int)*L*4,"cov_scales":None}
   inp=arm/"panel_stage2_input.json"; write(inp,{"n":len(y),"taus":TAUS,"q_tau":q,"baseline":base,"beta":beta.T,"se_tau":se.T,"n_eff":d,**cov}); outp=arm/"panel_stage2_output.json"; cmd=["bash",os.environ["MIXTURE_SCRIPT_PATH"],"--panel-stage2","--input",str(inp),"--output",str(outp),"--source",source,"--package-root",str(PINNED_PACKAGE_ROOT),"--package-commit",PINNED_PACKAGE_COMMIT,"--r-lib",env["r_lib"]]
   with (arm/"panel-stage2.stdout").open("w") as fo,(arm/"panel-stage2.stderr").open("w") as fe: subprocess.run(cmd,check=True,timeout=1100.,stdout=fo,stderr=fe)
   rr=json.loads(outp.read_text()); names=rr["parameter_names"]; rows=[]
   for one in rr["rows"]:
    got={canon(n):(float(e),float(s)) for n,e,s in zip(names,one["estimates"],one["ses"])}
    for p in PARAMS: rows.append({"snp_index":int(one["snp_index"]),"parameter":p,"estimate":got[p][0],"se":got[p][1],"status":"ok","reason":None})
   out={"status":"complete","attempted_snp_fits":M,"successful_snp_fits":M,"failures":0,"fit_seconds":time.perf_counter()-started,"rows":rows,"cost_label":"G-MoQ panel Stage 1 once plus pinned public param_gwas batch"}
  print(json.dumps(safe(out),allow_nan=False)); return 0
 except Exception as e:
  out=bad(f"{type(e).__name__}: {e}"); print(json.dumps(safe(out),allow_nan=False)); traceback.print_exc(); return 1

def panel_arm(method,G,y,seed,n,M,root):
 base=root/"panels"/f"n{n}_m{M}_rep{seed}"; base.mkdir(parents=True,exist_ok=True); payload=base/"input.pkl"
 if not payload.exists(): payload.write_bytes(pickle.dumps({"G":G,"y":y},protocol=pickle.HIGHEST_PROTOCOL))
 d=root/"panel_runs"/f"n{n}_m{M}_rep{seed}"/method; d.mkdir(parents=True,exist_ok=False); tf=d/"time-v.txt"; so=d/"stdout.log"; se=d/"stderr.log"; cmd=["/usr/bin/time","-v","-o",str(tf),"bash",os.environ["MIXTURE_SCRIPT_PATH"],"--panel-child",method,"--workdir",str(d),"--seed",str(seed),"--panel-input",str(payload)]
 try:
  with so.open("w") as fo,se.open("w") as fe:
   q=subprocess.Popen(cmd,stdout=fo,stderr=fe,start_new_session=True)
   try: q.communicate(timeout=7200.); rc=q.returncode; status="success" if rc==0 else "failed"
   except subprocess.TimeoutExpired: os.killpg(q.pid,signal.SIGKILL); q.wait(); rc=None; status="timeout"
 except Exception: rc=None; status="launch_failure"
 raw={}
 try: lines=so.read_text().splitlines(); raw=json.loads(lines[-1]) if lines else {}
 except Exception: pass
 return raw,dict(timeinfo(tf),process_status=status,return_code=rc,arm_status=status,log_dir=str(d))

def append_progress(path,phase,ident,n,M,seed,status,attempted,successful,failures):
 fields=["phase","id","n","M","seed","status","attempted_snp_fits","successful_snp_fits","failures"]; new=not Path(path).exists()
 with Path(path).open("a",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields); w.writeheader() if new else None; w.writerow(dict(zip(fields,[phase,ident,n,M,seed,status,attempted,successful,failures])))

def row(base,method,raw,ti,pop,case,effect,status_override=None,reason_override=None):
 truth_ok=status_override not in ("truth_failure",); target={canon(n):float(v) for n,v in zip(pop["theta_names"],pop["theta"])} if truth_ok else {p:None for p in PARAMS}; gen={"gamma":float(case["gamma"]),"delta1":float(case["delta1"]),"delta2":float(case["delta2"])}; out=[]
 fallback=status_override or ti.get("arm_status") or "failed"
 for rr in (raw.get("rows") or [{"parameter":p,"estimate":None,"se":None,"status":fallback,"reason":reason_override} for p in PARAMS]):
  if ti.get("process_status") not in (None,"success"): rr={**rr,"status":fallback}
  if rr.get("status")=="ok" and (rr.get("estimate") is None or rr.get("se") is None or not np.isfinite(rr["estimate"]) or not np.isfinite(rr["se"]) or rr["se"]<=0): rr={**rr,"status":"invalid_estimate_or_se"}
  p=canon(rr["parameter"]); good=rr.get("status")=="ok" and rr.get("estimate") is not None and rr.get("se") is not None and np.isfinite(rr["estimate"]) and np.isfinite(rr["se"])
  out.append({**base,"method":method,"parameter":p,"estimate":float(rr["estimate"]) if good else None,"se":float(rr["se"]) if good else None,"known_truth":target[p] if method.startswith("gmoq") else (gen[p] if truth_ok else None),"gmoq_projection_target":target[p],"ml_scientific_target":gen[p],"generating_effect":gen[p],"status":rr.get("status",fallback) if status_override is None else status_override,"reason":rr.get("reason") or reason_override,"fit_seconds":raw.get("fit_seconds"),"wall_seconds":ti.get("wall_seconds"),"peak_rss_kib":ti.get("peak_rss_kib"),"process_status":ti.get("process_status"),"return_code":ti.get("return_code"),"cost_label":raw.get("cost_label")})
 return out

def main():
 ap=argparse.ArgumentParser(); ap.add_argument("run_id",nargs="?",default="mixture-comparison-20260917-01"); ap.add_argument("reps",nargs="?",type=int,default=20); ap.add_argument("n_list",nargs="?",default="10000,50000,400000"); ap.add_argument("--child",choices=METHODS); ap.add_argument("--panel-child",choices=METHODS); ap.add_argument("--panel-size",type=int,default=100); ap.add_argument("--panel-reps",type=int,default=3); ap.add_argument("--workdir"); ap.add_argument("--input"); ap.add_argument("--panel-input"); ap.add_argument("--seed",type=int); a=ap.parse_args()
 if a.child: return child(argparse.Namespace(method=a.child,workdir=a.workdir,seed=a.seed,input=a.input))
 if a.panel_child: return panel_child(argparse.Namespace(method=a.panel_child,workdir=a.workdir,seed=a.seed,input=a.panel_input))
 ns=[int(x) for x in a.n_list.split(",") if x]; allowed={10000,50000,400000}
 if a.reps<1 or a.panel_reps<1 or a.panel_size<0 or not ns or any(n not in allowed for n in ns) or ns!=sorted(set(ns)): raise ValueError("reps positive, panel-size nonnegative, and N list must be increasing members of 10000,50000,400000")
 out=ROOT/".local"/a.run_id; out.mkdir(parents=True,exist_ok=False); rows=[]; failures=[]; panel_rows=[]; panel_timings=[]
 package_head=git(PINNED_PACKAGE_ROOT)
 if package_head!=PINNED_PACKAGE_COMMIT: raise RuntimeError(f"package is not pinned: {package_head}")
 prov={"schema_version":1,"run_id":a.run_id,"replicates":a.reps,"n_list":ns,"methods":list(METHODS),"panel_size":a.panel_size,"panel_reps":a.panel_reps,"panel_timeout_seconds":7200,"panel_methods":list(METHODS) if a.panel_size else [],"timeout_seconds":1200,"package_root":str(PINNED_PACKAGE_ROOT),"package_commit":PINNED_PACKAGE_COMMIT,"paper_commit":git(ROOT),"environment_manifest":str(ENV),"environment_sha256":sha(ENV),"script_sha256":sha(os.environ["MIXTURE_SCRIPT_PATH"]),"host":os.uname().nodename,"slurm_job":os.environ.get("SLURM_JOB_ID"),"threads":1,"truth_note":"population projection computed before timed arms; generating effects and truth are never supplied to fit_mixture_ml","timing_note":"Primary S5 is matched-panel timing from panel_timings.csv; recovery-arm independent-dataset timing is S5_timing_independent sensitivity. GNU time -v fresh process; G-MoQ panel uses public Stage 1 once plus pinned batch param_gwas/child R, ML fits each SNP sequentially. Failed panels remain recorded and are excluded from the primary plot; timing is descriptive, not a universal speed claim"}; write(out/"provenance.json",prov)
 for n in ns:
  for i in range(a.reps):
   seed=96000000+i; effect=float(np.random.default_rng(95000000+i).uniform(-.1,.1)); case={"id":"mixture_left","family":"mixture","maf":.3,"pi0":.5,"mu1":-2.,"mu2":2.,"sd1":1.,"sd2":1.,"gamma":0.,"delta1":effect,"delta2":0.}; base={"run_id":a.run_id,"dataset_id":f"n{n}_rep{i+1:03d}","n":n,"replicate":i+1,"seed":seed,"effect_seed":95000000+i,"case_id":"mixture_left","effect":effect}
   try: pop=truth.population_truth(case,TAUS)
   except Exception as e:
    for m in METHODS: rs=row(base,m,{"rows":[]},{}, {"theta":[0,0,0],"theta_names":PARAMS},case,"","truth_failure",f"{type(e).__name__}: {e}"); rows.extend(rs); failures.extend(rs)
    append_progress(out/"progress.csv","recovery",base["dataset_id"],n,0,base["seed"],"truth_failure",3,0,3); continue
   g_y=truth.draw(case,n,seed); g,y=g_y["g"],g_y["y"]
   for m in METHODS:
    try: raw,ti=arm(m,g,y,pop,seed,out); rs=row(base,m,raw,ti,pop,case,effect); rows.extend(rs); failures.extend([r for r in rs if r["status"]!="ok"])
    except Exception as e:
     rs=row(base,m,{"rows":[]},{"process_status":"orchestrator_failure","arm_status":"orchestrator_failure"},pop,case,effect,"orchestrator_failure",str(e)); rows.extend(rs); failures.extend(rs)
   append_progress(out/"progress.csv","recovery",base["dataset_id"],n,0,base["seed"],"complete" if all(r["status"]=="ok" for r in rows[-9:]) else "failed",3,sum(r["status"]=="ok" for r in rows[-9:])//3,3-sum(r["status"]=="ok" for r in rows[-9:])//3)
   write(out/"progress.json",{"status":"running","n_completed":sum(r["n"]==n and r["replicate"]<=i+1 for r in rows)//9,"n_rows":len(rows),"last_n":n,"last_replicate":i+1})
 if a.panel_size:
  for n in ns:
   for i in range(a.panel_reps):
    M=a.panel_size; seed=97000000+i; panel_id=f"n{n}_m{M}_rep{i+1:03d}"; null_case={"family":"mixture","maf":.3,"pi0":.5,"mu1":-2.,"mu2":2.,"sd1":1.,"sd2":1.,"gamma":0.,"delta1":0.,"delta2":0.}
    try:
     d0=truth.draw(null_case,n,seed); G=np.random.default_rng(seed+7919).binomial(2,.3,size=(n,M)).astype(np.int8); y0=d0["y"]
     for m in METHODS:
      try: raw,ti=panel_arm(m,G,y0,seed,n,M,out)
      except Exception as e: raw={"status":"orchestrator_failure","reason":str(e),"attempted_snp_fits":M,"successful_snp_fits":0,"failures":M,"rows":[]}; ti={"process_status":"orchestrator_failure","arm_status":"orchestrator_failure"}
      attempted=int(raw.get("attempted_snp_fits",M)); successful=int(raw.get("successful_snp_fits",0)); ps=raw.get("status",ti.get("arm_status","failed")); complete=ps=="complete" and successful==M and ti.get("process_status")=="success"
      if ps=="complete" and not complete: ps="process_failure"
      panel_timings.append({"run_id":a.run_id,"panel_id":panel_id,"n":n,"M":M,"seed":seed,"method":m,"fit_seconds":raw.get("fit_seconds"),"wall_seconds":ti.get("wall_seconds"),"peak_rss_kib":ti.get("peak_rss_kib"),"attempted_snp_fits":attempted,"successful_snp_fits":successful,"failures":max(attempted-successful,int(raw.get("failures",0))),"panel_status":"complete" if complete else ps,"process_status":ti.get("process_status"),"return_code":ti.get("return_code"),"reason":raw.get("reason")})
      fallback=ps
      rrs=raw.get("rows") or [{"snp_index":j,"parameter":p,"estimate":None,"se":None,"status":fallback,"reason":raw.get("reason")} for j in range(M) for p in PARAMS]
      for rr in rrs:
       good=rr.get("status")=="ok" and rr.get("estimate") is not None and rr.get("se") is not None and np.isfinite(rr["estimate"]) and np.isfinite(rr["se"]) and rr["se"]>0
       panel_rows.append({"run_id":a.run_id,"panel_id":panel_id,"n":n,"M":M,"seed":seed,"method":m,"snp_index":rr["snp_index"],"parameter":rr["parameter"],"estimate":float(rr["estimate"]) if good else None,"se":float(rr["se"]) if good else None,"status":"ok" if good else rr.get("status",fallback),"reason":rr.get("reason")})
     append_progress(out/"progress.csv","panel",panel_id,n,M,seed,"complete" if all(r["panel_status"]=="complete" for r in panel_timings[-3:]) else "failed",3*M, sum(r["successful_snp_fits"] for r in panel_timings[-3:]),sum(r["failures"] for r in panel_timings[-3:]))
     write(out/"progress.json",{"status":"running","phase":"panel","last_panel":panel_id,"n_rows":len(rows),"n_panel_rows":len(panel_rows)})
    except Exception as e:
     append_progress(out/"progress.csv","panel",panel_id,n,M,seed,"generation_failure",M,0,M); write(out/f"panel_generation_failure-{panel_id}.json",{"panel_id":panel_id,"error":f"{type(e).__name__}: {e}"})
 fields=["run_id","dataset_id","n","replicate","seed","effect_seed","case_id","effect","method","parameter","estimate","se","known_truth","gmoq_projection_target","ml_scientific_target","generating_effect","status","reason","fit_seconds","wall_seconds","peak_rss_kib","process_status","return_code","cost_label"]
 with (out/"results.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
 with (out/"failures.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows([r for r in rows if r["status"]!="ok"])
 pt_fields=["run_id","panel_id","n","M","seed","method","fit_seconds","wall_seconds","peak_rss_kib","attempted_snp_fits","successful_snp_fits","failures","panel_status","process_status","return_code","reason"]
 with (out/"panel_timings.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=pt_fields); w.writeheader(); w.writerows(panel_timings)
 pe_fields=["run_id","panel_id","n","M","seed","method","snp_index","parameter","estimate","se","status","reason"]
 with (out/"panel_estimates.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=pe_fields); w.writeheader(); w.writerows(panel_rows)
 sums=[]
 for n in ns:
  for m in METHODS:
   for p in PARAMS:
    z=[r for r in rows if r["n"]==n and r["method"]==m and r["parameter"]==p]; good=[r for r in z if r["status"]=="ok" and r["estimate"] is not None and r["se"] is not None]; err=np.array([r["estimate"]-r["known_truth"] for r in good]); sums.append({"n":n,"method":m,"parameter":p,"denominator_attempted":len(z),"denominator_success":len(good),"failures":len(z)-len(good),"timeouts":sum(r["status"]=="timeout" for r in z),"fit_failures":sum(r["status"] in ("fit_failure","hessian_failure","boundary_failure") for r in z),"truth_failures":sum(r["status"]=="truth_failure" for r in z),"bias":float(np.mean(err)) if len(err) else None,"rmse":float(np.sqrt(np.mean(err**2))) if len(err) else None,"mean_se":float(np.mean([r["se"] for r in good])) if good else None})
 with (out/"summary.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=list(sums[0])); w.writeheader(); w.writerows(sums)
 write(out/"COMPLETE.json",{"status":"data_complete","rows":len(rows),"failure_rows":len(failures),"panel_timing_rows":len(panel_timings),"panel_estimate_rows":len(panel_rows),"note":"Known-truth recovery is primary. Panel nulls are timing/finite-SE checks only, not calibration evidence. Agreement is not used as validation; estimates without finite SE and failed panels are excluded from figures."})
 try: subprocess.run(["bash",os.environ["MIXTURE_SCRIPT_PATH"],"--plot",str(out/"results.csv"),str(out/"figures"),str(out/"panel_timings.csv")],check=True)
 except Exception as e: write(out/"plot_failure.json",{"error":str(e)}); raise
 return 0
if __name__=="__main__": raise SystemExit(main())
PY
