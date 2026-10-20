#!/usr/bin/env bash
# Pooled against jackknife standard errors at common and rare variants
#
# What: normal phenotype, no genotype effect or a log-variance (dispersion) effect c = 0.1,
#   at every combination of N = 10,000, 50,000, 400,000 and MAF 0.05, 0.01, 0.001 (18 cells),
#   200 data sets per cell. Each data set is fitted by Stage 1 with both the pooled covariance
#   and the 1,000-block jackknife (identical point estimates, checked), then by Stage 2 for
#   location and dispersion with each covariance. Records the Stage 1 and Stage 2 SE ratios,
#   bias against the exact population projection target, the empirical SD, coverage of 95%
#   intervals, and the realised allele and tail-carrier counts.
# Paper: Supplementary Fig. S6; Supplementary Table S2 (study 6)
# Input: none (simulated data); ../run_examples.py, ../truth.py, ../stage2.R; the G-MoQ
#   packages through the environment made by ../bootstrap.sh
# Output: <ROOT>/.local/<run id>/{summary.csv,stage1_summary.csv,stage1_tau.csv,
#   stage2_parameters.csv,manifest.json} and two working plots
# Runs: on a cluster, one core, about 2 hours:
#   sbatch validation/figures/rare_variant_se.sh rare-se-20260918-01 200
#   Seed of replicate r (from 0) in cell c (from 0, ordered N, then MAF, then scenario):
#   98,100,000 + 10,000 c + r. The #SBATCH --chdir line holds the authors' cluster path;
#   edit it, and the user settings below, for your cluster.
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=8G
#SBATCH --chdir=/user/work/fh6520/Gmoq_paper-figures-20260917
#SBATCH --output=.local/rare-se-%j.out

set -euo pipefail
## User settings
# Directory that contains validation/ (here: analysis/08_simulations/code).
ROOT=${GMOQ_SIM_ROOT:-/user/work/fh6520/Gmoq_paper-figures-20260917}
# Environment made by validation/bootstrap.sh.
BOOT=${GMOQ_SIM_ENV:-/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-bootstrap-20260917-02}
source "$BOOT/activate.sh"
: "${GMOQ_STAGE1_PYTHON:?bootstrap did not set GMOQ_STAGE1_PYTHON}"
: "${GMOQ_RSCRIPT_BIN:?bootstrap did not set GMOQ_RSCRIPT_BIN}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

RUN_ID=${1:-rare-se-20260918-01}
REPS=${2:-200}
OUT="$ROOT/.local/$RUN_ID"
"$GMOQ_STAGE1_PYTHON" - "$RUN_ID" "$REPS" "$OUT" "$ROOT" "$BOOT/environment.json" "$GMOQ_RSCRIPT_BIN" <<'PY'
from __future__ import annotations
import csv, json, os, re, sys, traceback, time
from pathlib import Path
import numpy as np

run_id, reps_s, out_s, root_s, env_s, rscript = sys.argv[1:]
if not re.fullmatch(r"[A-Za-z0-9_.-]+", run_id): raise ValueError("unsafe run id")
reps = int(reps_s)
if reps < 1: raise ValueError("replicates must be positive")
ROOT, OUT, ENV = Path(root_s), Path(out_s), Path(env_s)
if OUT.exists(): raise RuntimeError(f"exclusive output already exists: {OUT}")
OUT.mkdir(parents=True); (OUT / "datasets").mkdir()

def atomic_json(path, value):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=lambda x: np.asarray(x).tolist()) + "\n")
    os.replace(tmp, path)

def finite(x, name):
    x = float(x)
    if not np.isfinite(x): raise RuntimeError(f"non-finite {name}")
    return x

manifest = {"status":"running", "run_id":run_id, "replicates":reps,
            "n_cells":18, "datasets_planned":18*reps, "package_commit":"9c47e7ba56370e708810fff8c2c9ad3b723cbccf",
            "environment_manifest":str(ENV), "threads":1, "n_blocks":1000, "taus":45,
            "grid":{"n":[10000,50000,400000], "maf":[.05,.01,.001], "scenarios":["normal_null","normal_dispersion"]}}
atomic_json(OUT / "manifest.json", manifest)
try:
    sys.path.insert(0, str(ROOT))
    from validation import run_examples as ex
    environment = ex.load_environment_manifest(ENV, ex.PINNED_PACKAGE_ROOT)
    if Path(rscript).resolve() != Path(environment["rscript"]).resolve():
        raise RuntimeError("GMOQ_RSCRIPT_BIN does not match the bootstrap manifest")
    core, package_provenance, compute_sigma0 = ex.load_pinned_core(ex.PINNED_PACKAGE_ROOT, environment)
    truth = ex._load_truth_module(ROOT)
except Exception as exc:
    atomic_json(OUT / "setup_failure.json", {"error":str(exc), "traceback":traceback.format_exc()})
    raise

manifest["source_commit"] = ex._git_head(ROOT)
manifest["ratio_definition"] = "jackknife_SE_divided_by_pooled_SE"
atomic_json(OUT/"manifest.json",manifest)
T = np.linspace(.02, .98, 45)
EXPECTED_PARAMS = ("beta_mu", "beta_sigma2")
config_base = {"density_floor":1e-8, "density_method":"silverman_gaussian", "min_mac":5.0,
               "min_var":1e-8, "n_blocks":1000, "threads":1, "panel_size":16}
scenarios = {"normal_null":{"a":0.0,"c":0.0}, "normal_dispersion":{"a":0.0,"c":0.1}}
cells = []
for n in (10000, 50000, 400000):
    for maf in (.05, .01, .001):
        for scenario, effects in scenarios.items():
            case = {"id":scenario, "family":"normal", "maf":maf, "a":effects["a"], "c":effects["c"],
                    "gamma":0.0,"delta1":0.0,"delta2":0.0,"pi0":.5,"mu1":-2.0,"mu2":2.0,"sd1":1.0,"sd2":1.0}
            cells.append({"cell_index":len(cells), "n":n, "maf":maf, "scenario":scenario, "case":case, "generating":effects})
cache = {}
for cell in cells:
    key = (cell["maf"], cell["scenario"])
    if key not in cache: cache[key] = truth.population_truth(dict(cell["case"]), T)

stage1_path, stage2_path = OUT/"stage1_tau.csv", OUT/"stage2_parameters.csv"
with stage1_path.open("w", newline="") as f1, stage2_path.open("w", newline="") as f2:
    w1, w2 = csv.writer(f1), csv.writer(f2)
    w1.writerow(["cell_index","replicate","n","maf","scenario","tau","beta_pooled","beta_jackknife","se_pooled","se_jackknife","se_ratio","actual_mac","carrier_count","homozygote_count","carrier_below_q02","carrier_above_q98","phenotype_q02","phenotype_q98"])
    w2.writerow(["cell_index","replicate","n","maf","scenario","param","source","estimate","se","ci_lower","ci_upper","p_value","projection_target","generating_parameter","generating_target","projection_bias","actual_mac","paired_se_ratio","paired_point_diff"])
    stage1_ratios, stage2_stats = {}, {}
    total = ok = failed = 0; started = time.monotonic()

    def add_stat(key, estimate, se, bias, covered, mac, ratio):
        z = stage2_stats.setdefault(key, {"estimate":[],"se":[],"bias":[],"covered":[],"mac":[],"ratio":[]})
        z["estimate"].append(estimate); z["se"].append(se); z["bias"].append(bias); z["covered"].append(bool(covered)); z["mac"].append(mac)
        if ratio is not None: z["ratio"].append(ratio)

    for rep in range(reps):
        for cell in cells:
            total += 1
            ci, n, maf, scenario = cell["cell_index"], cell["n"], cell["maf"], cell["scenario"]
            ddir = OUT/"datasets"/f"cell-{ci:02d}"/f"rep-{rep:03d}"; ddir.mkdir(parents=True)
            record = {"cell_index":ci,"replicate":rep,"n":n,"maf":maf,"scenario":scenario,"seed":98100000+ci*10000+rep,"status":"running"}
            captured = {}; original_draw = truth.draw
            def capture_draw(case, nn, seed):
                drawn = original_draw(case, nn, seed); g, y = np.asarray(drawn["g"]), np.asarray(drawn["y"],float)
                q02, q98 = np.quantile(y, [.02,.98], method="median_unbiased"); carriers = g > 0
                captured.update({"actual_mac":int(g.sum()),"carrier_count":int(carriers.sum()),"homozygote_count":int((g==2).sum()),"carrier_below_q02":int(np.sum(carriers&(y<q02))),"carrier_above_q98":int(np.sum(carriers&(y>q98))),"phenotype_q02":float(q02),"phenotype_q98":float(q98)})
                return drawn
            try:
                truth.draw = capture_draw
                payload = ex._stage1_example(core, truth, cell["case"], n, record["seed"], T, dict(config_base,maf=maf), compute_sigma0, population=cache[(maf,scenario)])
                truth.draw = original_draw
                if len(captured) != 7: raise RuntimeError("truth.draw wrapper did not capture exactly one dataset")
                bp, bj = np.asarray(payload["beta_pooled"],float), np.asarray(payload["beta_jackknife"],float)
                sp, sj = np.asarray(payload["se_tau_pooled"],float), np.asarray(payload["se_tau_jackknife"],float)
                diff = float(np.max(np.abs(bp-bj)))
                if diff > 1e-10+1e-10*max(1.,float(np.max(np.abs(bp)))): raise RuntimeError(f"Stage 1 paired point estimates disagree: {diff:g}")
                ratio = sj/sp
                if not np.all(np.isfinite(ratio)) or np.any(ratio<=0): raise RuntimeError("invalid Stage 1 SE ratio")
                for j, tau in enumerate(T):
                    w1.writerow([ci,rep,n,maf,scenario,float(tau),bp[j],bj[j],sp[j],sj[j],ratio[j],*(captured[k] for k in ("actual_mac","carrier_count","homozygote_count","carrier_below_q02","carrier_above_q98","phenotype_q02","phenotype_q98"))])
                    stage1_ratios.setdefault((ci,j),[]).append(float(ratio[j]))
                f1.flush(); input_path = ddir/"stage1_summary.json"; payload["package"] = package_provenance; ex._write_json(input_path,payload)
                results, errors = {}, {}
                for source in ("pooled","jackknife"):
                    try:
                        results[source] = ex._run_stage2(rscript,ROOT/"validation/stage2.R",input_path,ddir/f"stage2-{source}.json",ddir,source,ex.PINNED_PACKAGE_ROOT,ex.PINNED_PACKAGE_COMMIT,ex.PINNED_PACKAGE_VERSION,ENV,environment["r_lib"],1200.0)
                    except Exception as exc: errors[source] = str(exc)
                ratios, diffs = {}, {}
                if len(results) == 2:
                    names = [tuple(results[s]["estimate"]["parameter_names"]) for s in ("pooled","jackknife")]
                    if names[0] != names[1] or names[0] != EXPECTED_PARAMS: raise RuntimeError("unexpected or mismatched Stage 2 parameter names")
                    diffs = {p:abs(float(results["pooled"]["estimate"]["estimate"][p])-float(results["jackknife"]["estimate"]["estimate"][p])) for p in names[0]}
                    if max(diffs.values()) > 1e-10+1e-10*max(1.,max(diffs.values())): raise RuntimeError(f"Stage 2 paired point estimates disagree: {max(diffs.values()):g}")
                    ratios = {p:finite(results["jackknife"]["estimate"]["se"][p],"SE")/finite(results["pooled"]["estimate"]["se"][p],"SE") for p in names[0]}
                    if any(v<=0 for v in ratios.values()): raise RuntimeError("non-positive paired Stage 2 SE ratio")
                rows = []
                for source, result in results.items():
                    est, target = result["estimate"], result["target"]["theta_pop"]
                    names = tuple(est["parameter_names"])
                    if names != EXPECTED_PARAMS: raise RuntimeError("unexpected Stage 2 parameter names")
                    for j, p in enumerate(names):
                        e, se = finite(est["estimate"][p],"estimate"), finite(est["se"][p],"SE")
                        if se <= 0: raise RuntimeError("Stage 2 SE is not positive")
                        lo, hi = finite(est["ci_lower"][p],"CI lower"), finite(est["ci_upper"][p],"CI upper")
                        pv = est.get("p_value",{}); pv = None if pv is None else pv.get(p); pv = "" if pv is None else finite(pv,"p-value")
                        proj, gen_name = finite(target[p],"projection target"), ("a","c")[j]; gen = float(cell["generating"][gen_name]); bias = e-proj
                        row = [ci,rep,n,maf,scenario,p,source,e,se,lo,hi,pv,proj,gen_name,gen,bias,captured["actual_mac"],ratios.get(p,""),diffs.get(p,"")]
                        rows.append((row,(n,maf,scenario,p,source),e,se,bias,lo<=proj<=hi,ratios.get(p),captured["actual_mac"]))
                for row,key,e,se,bias,cov,rr,mac in rows: w2.writerow(row); add_stat(key,e,se,bias,cov,mac,rr)
                f2.flush()
                record.update({"status":"ok" if not errors else "failed_stage2","actual":captured,"stage2_errors":errors,"stage2_se_ratios":ratios})
                if errors: ex._write_json(ddir/"failure.json",record); failed += 1
                else: ok += 1
            except Exception as exc:
                truth.draw = original_draw; failed += 1
                record.update({"status":"failed","error":str(exc),"traceback":traceback.format_exc(),"actual":captured})
                ex._write_json(ddir/"failure.json",record)
            ex._write_json(ddir/"record.json",record)
            atomic_json(OUT/"progress.json",{"status":"running","attempted":total,"success":ok,"failed":failed,"last":record,"elapsed_seconds":time.monotonic()-started})

    def summary_row(key, z, attempted):
        vals, ses, bias = map(np.asarray,(z["estimate"],z["se"],z["bias"])); ns=len(vals)
        emp=float(np.std(vals,ddof=1)) if ns>1 else None; meanse=float(np.mean(ses)) if ns else None; p=float(np.mean(z["covered"])) if ns else None
        return [*key,attempted,ns,attempted-ns,meanse,emp,"" if not emp else meanse/emp,float(np.mean(bias)) if ns else None,p,"" if not ns else np.sqrt(p*(1-p)/ns),float(np.mean(z["mac"])) if ns else None,float(np.median(z["ratio"])) if z["ratio"] else None,float(np.quantile(z["ratio"],.05)) if z["ratio"] else None,float(np.quantile(z["ratio"],.95)) if z["ratio"] else None]
    with (OUT/"summary.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["n","maf","scenario","param","source","attempted","success","fail","meanSE","empiricalSD","SE_over_empSD","projectionbias","coverage","coverageMCSE","actualMACmean","medianpairedSEratio","pairedSEratio_p05","pairedSEratio_p95"])
        for cell in cells:
            for p in EXPECTED_PARAMS:
                for source in ("pooled","jackknife"):
                    key=(cell["n"],cell["maf"],cell["scenario"],p,source); w.writerow(summary_row(key,stage2_stats.get(key,{"estimate":[],"se":[],"bias":[],"covered":[],"mac":[],"ratio":[]}),reps))
    with (OUT/"stage1_summary.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["n","maf","scenario","tau","attempted","success","fail","medianratio","ratio_p05","ratio_p95"])
        for cell in cells:
            for j,tau in enumerate(T):
                vals=stage1_ratios.get((cell["cell_index"],j),[]); w.writerow([cell["n"],cell["maf"],cell["scenario"],float(tau),reps,len(vals),reps-len(vals),float(np.median(vals)) if vals else None,float(np.quantile(vals,.05)) if vals else None,float(np.quantile(vals,.95)) if vals else None])
    manifest.update({"status":"ok" if failed==0 else "completed_with_failures","datasets_attempted":total,"datasets_success":ok,"datasets_failed":failed,"package":package_provenance,"environment":{"manifest":str(ENV),"rscript":rscript,"r_library":environment["r_lib"]},"cells":[{"cell_index":c["cell_index"],"n":c["n"],"maf":c["maf"],"scenario":c["scenario"],"generating":c["generating"],"projection_target":cache[(c["maf"],c["scenario"])]["theta"]} for c in cells]})
    atomic_json(OUT/"manifest.json",manifest); atomic_json(OUT/"progress.json",manifest)
    if failed: raise SystemExit(1)
PY

if [[ "${GMOQ_MAKE_FIGURE:-1}" == 1 ]]; then
  "$GMOQ_RSCRIPT_BIN" --vanilla - "$OUT/summary.csv" "$OUT" <<'RS'
args <- commandArgs(trailingOnly=TRUE); suppressPackageStartupMessages(library(ggplot2))
d <- read.csv(args[1], check.names=FALSE); d <- d[d$source == "pooled" & d$success > 0, ]; d$cell <- interaction(d$n, d$scenario, sep=" / ")
p <- ggplot(d, aes(factor(maf), medianpairedSEratio, colour=param, group=param)) + geom_point() + geom_line() + geom_hline(yintercept=1, linetype=2) + facet_wrap(~cell, scales="free_x") + labs(x="MAF", y="Median jackknife / pooled SE") + theme_bw()
ggsave(file.path(args[2], "stage2_se_ratios.png"), p, width=11, height=7, dpi=160)
d <- read.csv(args[1], check.names=FALSE); d <- d[d$success > 0, ]; d$cell <- interaction(d$n, d$scenario, d$param, sep=" / ")
p <- ggplot(d, aes(factor(maf), coverage, colour=source, group=source)) + geom_point() + geom_line() + geom_hline(yintercept=.95, linetype=2) + facet_wrap(~cell, scales="free_x") + labs(x="MAF", y="95% CI coverage of projection target") + theme_bw()
ggsave(file.path(args[2], "stage2_coverage.png"), p, width=11, height=7, dpi=160)
RS
fi
