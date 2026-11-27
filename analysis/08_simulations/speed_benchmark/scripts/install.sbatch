#!/bin/bash
#SBATCH --job-name=sim-speed-build
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -eo pipefail
case "$(hostname)" in bp1-compute*) ;; *) exit 2 ;; esac
hostname
source /etc/profile.d/modules.sh
module load git/2.45.1-pqk5 languages/R/4.5.1
set -u
pin=166635d49527e27f8a714247d551c17093e02a55
root=/user/work/fh6520/gmoq-paper-v1/pkg/gmoq-166635d4
lib=/user/work/fh6520/gmoq-paper-v1/Rlib/gmoq-166635d4
export CONDA_PREFIX=/user/home/fh6520/miniforge3/envs/fungwas-stage1
export PATH="$CONDA_PREFIX/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py="$CONDA_PREFIX/bin/python"
cd "$root"
test "$(git rev-parse HEAD)" = "$pin"
test -z "$(git status --porcelain --untracked-files=no)"
mkdir -p "$lib"
(cd stage1 && "$py" setup.py build_ext --inplace)
export PYTHONPATH="$root/stage1"
"$py" -c 'from gmoq import core; assert core.cpp_ext.openmp_info()["max_threads"]==1; print(core.cpp_ext.__file__)'
R_LIBS_USER=/user/work/fh6520/r_packages R CMD INSTALL --no-test-load -l "$lib" stage2
export RERUN_R_LIB="$lib"
version=$(Rscript -e '.libPaths(c(Sys.getenv("RERUN_R_LIB"), "/user/work/fh6520/r_packages", .libPaths())); suppressPackageStartupMessages(library(gmoq)); cat(as.character(packageVersion("gmoq")))')
test "$version" = "0.8.1"
extension_hash=$(sha256sum stage1/gmoq/_stage1_cpp*.so | cut -c1-64)
cat > "$lib/PIN.json.tmp" <<JSON
{"commit":"$pin","stage1_root":"$root","stage1_python":"$py","stage1_extension_sha256":"$extension_hash","r_package":"gmoq","r_package_version":"$version","r_lib":"$lib","r_module":"languages/R/4.5.1","install_job":"$SLURM_JOB_ID"}
JSON
mv "$lib/PIN.json.tmp" "$lib/PIN.json"
cat "$lib/PIN.json"
