#!/usr/bin/env python3
"""Stage 1 against ordinary RIF-OLS: coefficients, standard errors and timing.

What: simulates four workloads (N = 2,000 with 1 and 256 variants; N = 400,000 with 1 and
  1,000 variants; phenotype 0.2 x centred genotype of the first variant plus N(0, 1) noise,
  genotypes binomial(2, 0.3)), then fits each workload with three methods in fresh,
  single-core processes, three timed repeats each (36 runs): G-MoQ Stage 1 with the pooled
  covariance, G-MoQ Stage 1 with the 1,000-block jackknife, and ordinary RIF-OLS (numpy
  lstsq on the same RIFs, one regression per quantile level). It then checks that the
  coefficients agree and exports the comparison and timing tables.
Paper: Supplementary Figs. S1, S2; Supplementary Table S2 (studies 1 and 2)
Input: none (simulated data); the G-MoQ packages, through the environment made by
  ../bootstrap.sh
Output: a new run directory under .local/ (beside validation/) holding tables/{S1_pairs,S1_SE,S1_checks,
  S2_timings}.csv, verification.json, and draft plots
Runs: on a cluster, one core, about 25 minutes, through stage1.sbatch:
  python validation/figures/run_stage1.py --attempt 05 --n-large 400000
  Seeds: data 940001-940004 (one per workload); run order 940100.
  The default environment path below is the authors' cluster path; set
  GMOQ_VALIDATION_ENV (exported by the environment's activate.sh) to override it.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from validation.benchmark import run_benchmark
from validation.run_examples import PINNED_PACKAGE_COMMIT
# User setting: the environment made by ../bootstrap.sh (its activate.sh exports
# GMOQ_VALIDATION_ENV). The default is the authors' cluster path.
ENV = Path(os.environ.get('GMOQ_VALIDATION_ENV',
                          '/user/work/fh6520/Gmoq_paper-validation-20260917/.local/validation-bootstrap-20260917-02')) / 'environment.json'
WORKLOADS = [(2000, 1, 940001), (2000, 256, 940002)]
LARGE_WORKLOAD_SEEDS = [(1, 940003), (1000, 940004)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--attempt', default='01')
    ap.add_argument('--n-large', type=int, choices=(50000, 400000), default=50000)
    args = ap.parse_args()
    if len(args.attempt) != 2 or not args.attempt.isdigit():
        raise ValueError('attempt must be two digits')
    experiment = 'stage1-figures-20260917-' + args.attempt
    out = ROOT / '.local' / experiment
    out.mkdir(parents=True, exist_ok=False)
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    os.environ['VALIDATION_COMPUTE_HOST'] = '1'
    topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
    provenance = {'paper_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'package_commit': PINNED_PACKAGE_COMMIT, 'host': socket.gethostname(),
                  'job': os.environ.get('SLURM_JOB_ID', ''), 'logical_cpu': cpu,
                  'core_id': (topology / 'core_id').read_text().strip(),
                  'socket_id': (topology / 'physical_package_id').read_text().strip(),
                  'platform': platform.platform(), 'environment_manifest': str(ENV),
                  'scope': 'Stage1 coefficient and SE calculation; pooled and jackknife separate; not conditional quantile regression'}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    methods_script = ROOT / 'validation/figures/stage1_methods.py'
    workloads = []
    workload_specs = WORKLOADS + [(args.n_large, m, seed) for m, seed in LARGE_WORKLOAD_SEEDS]
    for n, m, seed in workload_specs:
        wid = f'n{n}_m{m}'
        data = out / 'inputs' / wid
        subprocess.run([sys.executable, str(methods_script), 'generate', '--n', str(n), '--m', str(m),
                        '--seed', str(seed), '--output-dir', str(data)], check=True)
        workloads.append({'workload_id': wid, 'label': f'N={n}, M={m}',
                          'execution_mode': 'isolated_single_snp' if m == 1 else 'amortized_scan',
                          'placeholders': {'input_path': str(data)}, 'timing_file': 'timing.json'})
    methods = []
    for method in ("gmoq", "ordinary_rif_ols", "gmoq_jackknife"):
        methods.append({'method_id': method, 'label': method,
                        'command': [sys.executable, str(methods_script), 'fit', '--input-dir', '{input_path}',
                                    '--method', method, '--output-dir', '{output_dir}', '--threads', '{threads}',
                                    '--package-commit', '{package_pin}', '--environment-manifest', str(ENV)]})
    manifest = {'schema_version': 'validation-benchmark-manifest-1', 'experiment_id': experiment,
                'package_pin': PINNED_PACKAGE_COMMIT, 'output_root': str(out / 'benchmark'),
                'compute_host_only': True, 'physical_cores': 1, 'physical_cpu_ids': [cpu], 'threads': 1,
                'timeout_seconds': 1800 if args.n_large == 400000 else 600,
                'interleave_seed': 940100, 'warmup_count': 0, 'repeats': 3,
                'workloads': workloads, 'methods': methods}
    (out / 'command_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    results = run_benchmark(manifest)
    if len(results) != 36 or any(row['status'] != 'success' for row in results):
        raise RuntimeError('benchmark arm failed; retained results/logs; do not plot as success')
    subprocess.run([sys.executable, str(ROOT / 'validation/figures/export_stage1.py'),
                    '--benchmark-dir', str(out / 'benchmark'), '--out-dir', str(out / "tables"),
                    '--n-large', str(args.n_large)], check=True)
    subprocess.run([os.environ["GMOQ_RSCRIPT_BIN"], "--vanilla", str(ROOT / "validation/figures/plot_stage1.R"), str(out / "tables"), str(out / "figures")], check=True)
    (out / "COMPLETE").write_text('S1/S2 pipeline completed; inspect the figures.\n')
    print(json.dumps({'status': 'complete', 'output': str(out), 'arms': len(results)}))


if __name__ == '__main__':
    main()
