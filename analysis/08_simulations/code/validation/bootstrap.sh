#!/usr/bin/env bash
# Build the software environment used by every simulation study
#
# What: checks that a checkout of the G-MoQ package is at the pinned Stage 1 commit with no
#   local changes, compiles the Stage 1 extension in the gmoq-stage1 conda environment,
#   installs Stage 2 (gmoq) into a new private R library, and writes
#   environment.json (interpreter paths, versions and source hashes, read by
#   run_examples.load_environment_manifest) and activate.sh (exports GMOQ_STAGE1_PYTHON,
#   GMOQ_RSCRIPT_BIN and GMOQ_VALIDATION_ENV for the drivers). It installs no other
#   dependencies.
# Paper: Supplementary Figs. S1-S6; Supplementary Tables S1, S2 (setup for studies 1-7)
# Input: a G-MoQ package checkout; a conda environment gmoq-stage1 with numpy, scipy and a
#   C++ compiler; R 4.5 with jsonlite and the Stage 2 dependencies
# Output: <out>/{environment.json,activate.sh,Rlib/}
# Runs: once, inside a cluster compute-node job, a few minutes:
#   bash validation/bootstrap.sh --out <new directory>
#   The settings below, and the module and conda lines further down (git/2.45.1-pqk5,
#   languages/R/4.5.1, the conda.sh path), are the authors' cluster settings; edit them for
#   your system.
set -Eeuo pipefail

## User settings
readonly PACKAGE_ROOT="/user/work/fh6520/G-MoQ-validation-9c47e7b"
readonly PACKAGE_PIN="9c47e7ba56370e708810fff8c2c9ad3b723cbccf"
readonly CONDA_SH="/user/home/fh6520/miniforge3/etc/profile.d/conda.sh"
readonly GIT_MODULE="git/2.45.1-pqk5"
readonly R_MODULE="languages/R/4.5.1"
OUT=""

usage() {
    cat >&2 <<'EOF'
usage: bootstrap.sh --out /new/environment/output/path

The output path must not already exist and this script must run inside a Slurm
compute-node allocation.
EOF
}

die() {
    echo "bootstrap: ERROR: $*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --out)
            (($# >= 2)) || { usage; exit 2; }
            OUT=$2
            shift 2
            ;;
        --out=*)
            OUT=${1#--out=}
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            die "unknown argument: $1"
            ;;
    esac
done

[[ -n "$OUT" ]] || { usage; die "--out is required"; }
[[ "$OUT" != "/" && "$OUT" != "." && "$OUT" != ".." ]] || die "refusing unsafe output path: $OUT"

# Do not allow a login-node execution, even if a caller accidentally exports a
# misleading path.  SLURM_JOB_ID is required for an actual allocation.
[[ -n "${SLURM_JOB_ID:-}" ]] || die "no SLURM_JOB_ID; run this on a compute-node allocation, not a login node"
HOST=$(hostname -s)
case "$HOST" in
    login*|*login*|submit*|*submit*|head*|*gateway*)
        die "host $HOST looks like a login/submit host; refusing to build there"
        ;;
esac

OUT_PARENT=$(dirname -- "$OUT")
mkdir -p -- "$OUT_PARENT"
OUT_PARENT=$(cd -- "$OUT_PARENT" && pwd -P)
OUT="$OUT_PARENT/$(basename -- "$OUT")"
[[ ! -e "$OUT" && ! -L "$OUT" ]] || die "output path already exists (must be unique and new): $OUT"
mkdir -- "$OUT"

on_exit() {
    local rc=$?
    if ((rc != 0)); then
        printf 'bootstrap failed on %s with status %d\n' "$HOST" "$rc" > "$OUT/FAILED"
    fi
    exit "$rc"
}
trap on_exit EXIT

printf 'bootstrap: host=%s slurm_job=%s output=%s\n' "$HOST" "$SLURM_JOB_ID" "$OUT"
printf 'bootstrap: package=%s pin=%s\n' "$PACKAGE_ROOT" "$PACKAGE_PIN"

# Compute images do not necessarily put Git on PATH.  Initialise Environment
# Modules before any Git check, then load the pinned Git module on its own line.
if ! type module >/dev/null 2>&1; then
    for MODULE_INIT in /usr/share/Modules/init/bash /etc/profile.d/modules.sh /usr/local/Modules/init/bash; do
        if [[ -r "$MODULE_INIT" ]]; then
            # shellcheck disable=SC1090
            source "$MODULE_INIT"
            break
        fi
    done
fi
if ! type module >/dev/null 2>&1; then
    die "Environment Modules command is unavailable; cannot load $GIT_MODULE"
fi
module load git/2.45.1-pqk5

command -v git >/dev/null 2>&1 || die "git is required to verify the pinned source"
command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required to record source hashes"
[[ -d "$PACKAGE_ROOT" ]] || die "package source directory does not exist: $PACKAGE_ROOT"
[[ -f "$PACKAGE_ROOT/stage1/setup.py" ]] || die "missing Stage 1 setup.py under pinned source"
[[ -f "$PACKAGE_ROOT/stage2/DESCRIPTION" ]] || die "missing Stage 2 DESCRIPTION under pinned source"
[[ -f "$PACKAGE_ROOT/stage1/src/stage1_fast.hpp" ]] || die "pinned fast-kernel header is missing"
[[ -f "$PACKAGE_ROOT/stage1/src/stage1_fast_hot.inc" ]] || die "pinned fast-kernel include is missing"

HEAD=$(git -C "$PACKAGE_ROOT" rev-parse HEAD)
[[ "$HEAD" == "$PACKAGE_PIN" ]] || die "package HEAD is $HEAD, expected pinned SHA $PACKAGE_PIN"
TRACKED_STATUS=$(git -C "$PACKAGE_ROOT" status --porcelain=v1 --untracked-files=no)
[[ -z "$TRACKED_STATUS" ]] || die "package has tracked modifications; refusing to build:\n$TRACKED_STATUS"
# The archive hash covers the exact committed source, independently of the Git
# object name, and is recorded in the reusable environment manifest.
SOURCE_TREE_SHA=$(git -C "$PACKAGE_ROOT" rev-parse "${PACKAGE_PIN}^{tree}")
SOURCE_ARCHIVE_SHA=$(git -C "$PACKAGE_ROOT" archive --format=tar "$PACKAGE_PIN" | sha256sum | awk '{print $1}')
printf 'bootstrap: verified clean tracked source at %s (tree %s)\n' "$HEAD" "$SOURCE_TREE_SHA"

[[ -r "$CONDA_SH" ]] || die "conda initialisation script is missing: $CONDA_SH"
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
if ! type conda >/dev/null 2>&1; then
    die "conda command was not defined by $CONDA_SH"
fi
conda activate gmoq-stage1
[[ "${CONDA_DEFAULT_ENV:-}" == "gmoq-stage1" ]] || die "conda activation did not select gmoq-stage1"
PYTHON_BIN=$(type -P python || true)
[[ -n "$PYTHON_BIN" && -x "$PYTHON_BIN" ]] || die "python was not found after activating gmoq-stage1"
case "$PYTHON_BIN" in
    "${CONDA_PREFIX:-}"/*) ;;
    *) die "python is not from gmoq-stage1: $PYTHON_BIN" ;;
esac
export PYTHONPATH="$PACKAGE_ROOT/stage1${PYTHONPATH:+:$PYTHONPATH}"

# This is deliberately the pinned source build command; it may create only
# compiled/untracked artifacts in stage1, not package-source edits.
cd "$PACKAGE_ROOT/stage1"
python setup.py build_ext --inplace
cd "$PACKAGE_ROOT"

HEAD_AFTER=$(git -C "$PACKAGE_ROOT" rev-parse HEAD)
[[ "$HEAD_AFTER" == "$PACKAGE_PIN" ]] || die "package HEAD changed during build: $HEAD_AFTER"
TRACKED_STATUS_AFTER=$(git -C "$PACKAGE_ROOT" status --porcelain=v1 --untracked-files=no)
[[ -z "$TRACKED_STATUS_AFTER" ]] || die "Stage 1 build modified tracked package files:\n$TRACKED_STATUS_AFTER"

"$PYTHON_BIN" - <<'PY' > "$OUT/python-info.json"
import importlib
import json
import os
import sys

required = ("numpy", "scipy", "pytest")
records = {}
missing = []
for name in required:
    try:
        module = importlib.import_module(name)
        records[name] = {
            "version": str(getattr(module, "__version__", "unknown")),
            "path": os.path.realpath(getattr(module, "__file__", "")),
        }
    except Exception as exc:  # pragma: no cover - executed by bootstrap
        missing.append(f"{name}: {type(exc).__name__}: {exc}")
if missing:
    raise SystemExit(
        "Missing/unusable Python packages (no downloads attempted):\n  "
        + "\n  ".join(missing)
        + "\npython=" + sys.executable
        + "\nsys.path=" + repr(sys.path)
    )

try:
    package = importlib.import_module("gmoq")
    core = importlib.import_module("gmoq.core")
except Exception as exc:
    raise SystemExit(
        "Pinned Stage 1 package import failed (compiled extension/dependencies): "
        + f"{type(exc).__name__}: {exc}\npython={sys.executable}"
    )

if not bool(getattr(core, "HAVE_CPP", False)):
    raise SystemExit(
        "Stage 1 C++ extension did not load; refusing the NumPy fallback. "
        "Rebuild the pinned source on a node with its compiler/Armadillo toolchain."
    )
extension = getattr(core, "cpp_ext", None)
extension_path = os.path.realpath(getattr(extension, "__file__", ""))
if not extension_path or not os.path.isfile(extension_path):
    raise SystemExit(f"Stage 1 extension path is missing or not a file: {extension_path!r}")

records["gmoq"] = {
    "version": str(getattr(package, "__version__", "unknown")),
    "path": os.path.realpath(getattr(package, "__file__", "")),
}
records["gmoq.core"] = {
    "path": os.path.realpath(getattr(core, "__file__", "")),
    "have_cpp": True,
}
print(json.dumps({
    "executable": os.path.realpath(sys.executable),
    "conda_prefix": os.path.realpath(os.environ.get("CONDA_PREFIX", "")),
    "python_version": sys.version,
    "modules": records,
    "extension": {"path": extension_path, "module": "gmoq._stage1_cpp"},
}, indent=2, sort_keys=True))
PY

EXTENSION_PATH=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["extension"]["path"])' "$OUT/python-info.json")
[[ -f "$EXTENSION_PATH" ]] || die "compiled Stage 1 extension is missing: $EXTENSION_PATH"
EXTENSION_SHA=$(sha256sum -- "$EXTENSION_PATH" | awk '{print $1}')

# Do not let the conda environment's optional R installation leak into the
# module R ABI.  The module load below is intentionally a standalone command,
# not a subshell or a pipeline.  Modules were initialised and Git loaded before
# the source checks above.
unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE || true
module load languages/R/4.5.1
hash -r
R_BIN=$(type -P R || true)
RSCRIPT_BIN=$(type -P Rscript || true)
[[ -n "$R_BIN" && -x "$R_BIN" ]] || die "R was not found after loading $R_MODULE"
[[ -n "$RSCRIPT_BIN" && -x "$RSCRIPT_BIN" ]] || die "Rscript was not found after loading $R_MODULE"
R_BIN=$(readlink -f "$R_BIN")
RSCRIPT_BIN=$(readlink -f "$RSCRIPT_BIN")
case "$R_BIN" in
    "${CONDA_PREFIX:-}"/*) die "R resolves into conda ($R_BIN); refusing a mixed R ABI" ;;
esac
case "$RSCRIPT_BIN" in
    "${CONDA_PREFIX:-}"/*) die "Rscript resolves into conda ($RSCRIPT_BIN); refusing a mixed R ABI" ;;
esac

RLIB="$OUT/Rlib"
mkdir -- "$RLIB"
export R_LIBS_USER="$RLIB"

# Check dependencies before installing, so a missing module package fails with
# an explicit list rather than triggering any implicit installation behaviour.
"$RSCRIPT_BIN" --vanilla - "$RLIB" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
rlib <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(unique(c(rlib, .libPaths())))
required <- c("Rcpp", "RcppArmadillo", "data.table", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    cat("Missing required R packages (no downloads attempted): ",
        paste(missing, collapse = ", "), "\n", sep = "", file = stderr())
    cat("R library paths:\n  ", paste(.libPaths(), collapse = "\n  "), "\n", sep = "", file = stderr())
    quit(status = 1)
}
cat("R dependency preflight passed; optional dglm=",
    requireNamespace("dglm", quietly = TRUE), "\n", sep = "")
RS

# Install the Stage 2 package from the same pinned, clean source tree into the
# new private library.  No shared installation can satisfy this package check.
R_LIBS_USER="$RLIB" "$R_BIN" CMD INSTALL --no-multiarch --library="$RLIB" "$PACKAGE_ROOT/stage2"

"$RSCRIPT_BIN" --vanilla - "$RLIB" "$R_BIN" "$RSCRIPT_BIN" "$OUT/r-info.json" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
rlib <- normalizePath(args[[1]], mustWork = TRUE)
rbin <- args[[2]]
rsbin <- args[[3]]
out <- args[[4]]
.libPaths(unique(c(rlib, .libPaths())))
required <- c("Rcpp", "RcppArmadillo", "data.table", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    stop("Required R packages became unavailable: ", paste(missing, collapse = ", "))
}

stage2_path <- find.package("gmoq", lib.loc = rlib, quiet = TRUE)
if (!length(stage2_path)) {
    stop("Pinned gmoq was not found in isolated Rlib: ", rlib)
}
stage2_path <- normalizePath(stage2_path[[1]], mustWork = TRUE)
rlib_prefix <- paste0(rlib, .Platform$file.sep)
if (!startsWith(stage2_path, rlib_prefix)) {
    stop("gmoq resolved outside isolated Rlib: ", stage2_path)
}
suppressPackageStartupMessages(library(gmoq, lib.loc = rlib))
if (as.character(packageVersion("gmoq")) != "0.4.0") {
    stop("Pinned Stage 2 version is not 0.4.0: ", packageVersion("gmoq"))
}

package_info <- function(name) {
    path <- find.package(name, quiet = TRUE)
    list(
        version = as.character(packageVersion(name)),
        path = normalizePath(path[[1]], mustWork = TRUE)
    )
}
optional_dglm <- requireNamespace("dglm", quietly = TRUE)
info <- list(
    executables = list(R = rbin, Rscript = rsbin),
    version = R.version.string,
    platform = R.version$platform,
    R_home = R.home(),
    library_paths = .libPaths(),
    packages = list(
        Rcpp = package_info("Rcpp"),
        RcppArmadillo = package_info("RcppArmadillo"),
        data.table = package_info("data.table"),
        jsonlite = package_info("jsonlite"),
        dglm = if (optional_dglm) package_info("dglm") else list(available = FALSE)
    ),
    gmoq = list(
        version = as.character(packageVersion("gmoq")),
        installed_path = stage2_path,
        isolated_library = rlib
    )
)
jsonlite::write_json(info, out, auto_unbox = TRUE, pretty = TRUE)
RS

STAGE2_DESCRIPTION="$RLIB/gmoq/DESCRIPTION"
[[ -f "$STAGE2_DESCRIPTION" ]] || die "installed Stage 2 DESCRIPTION is missing from isolated Rlib"
STAGE2_DESCRIPTION_SHA=$(sha256sum -- "$STAGE2_DESCRIPTION" | awk '{print $1}')

# Later scripts source this file instead of relying on whichever PATH or R
# library happened to be present in their submission shell.
ACTIVATE="$OUT/activate.sh"
{
    cat <<'ACT'
#!/usr/bin/env bash
# Source on a compute node: source /path/to/activate.sh
if ! type module >/dev/null 2>&1; then
    for _module_init in /usr/share/Modules/init/bash /etc/profile.d/modules.sh /usr/local/Modules/init/bash; do
        if [[ -r "$_module_init" ]]; then
            # shellcheck disable=SC1090
            source "$_module_init"
            break
        fi
    done
fi
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate gmoq-stage1
module load git/2.45.1-pqk5
unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE || true
module load languages/R/4.5.1
ACT
    printf 'export GMOQ_VALIDATION_ENV=%q\n' "$OUT"
    printf 'export GMOQ_PACKAGE_ROOT=%q\n' "$PACKAGE_ROOT"
    printf 'export GMOQ_PACKAGE_PIN=%q\n' "$PACKAGE_PIN"
    printf 'export GMOQ_STAGE1_PYTHON=%q\n' "$PYTHON_BIN"
    printf 'export GMOQ_STAGE2_RLIB=%q\n' "$RLIB"
    printf 'export GMOQ_R_BIN=%q\n' "$R_BIN"
    printf 'export GMOQ_RSCRIPT_BIN=%q\n' "$RSCRIPT_BIN"
    printf 'export R_LIBS_USER=%q\n' "$RLIB"
    printf 'export PYTHONPATH=%q:${PYTHONPATH:-}\n' "$PACKAGE_ROOT/stage1"
} > "$ACTIVATE"
chmod 0555 "$ACTIVATE"

"$PYTHON_BIN" - "$OUT/python-info.json" "$OUT/r-info.json" "$OUT/environment.json" \
    "$PACKAGE_ROOT" "$PACKAGE_PIN" "$HEAD" "$SOURCE_TREE_SHA" "$SOURCE_ARCHIVE_SHA" \
    "$EXTENSION_PATH" "$EXTENSION_SHA" "$STAGE2_DESCRIPTION_SHA" "$ACTIVATE" "$RLIB" <<'PY'
import json
import sys
from datetime import datetime, timezone

(
    python_info_path, r_info_path, manifest_path,
    package_root, package_pin, head, tree_sha, archive_sha,
    extension_path, extension_sha, stage2_description_sha, activate, rlib,
) = sys.argv[1:]
with open(python_info_path, encoding="utf-8") as handle:
    python_info = json.load(handle)
with open(r_info_path, encoding="utf-8") as handle:
    r_info = json.load(handle)
manifest = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "requested_environment": {
        "conda": "gmoq-stage1",
        "git_module": "git/2.45.1-pqk5",
        "r_module": "languages/R/4.5.1",
    },
    "package": {
        "root": package_root,
        "commit": package_pin,
        "head_verified": head,
        "tree_sha": tree_sha,
        "git_archive_sha256": archive_sha,
        "tracked_clean_before_build": True,
        "tracked_clean_after_build": True,
        "stage1_fast_kernel_files": [
            "stage1/src/stage1_fast.hpp",
            "stage1/src/stage1_fast_hot.inc",
        ],
    },
    "python": python_info,
    "stage1_extension": {
        "path": extension_path,
        "sha256": extension_sha,
    },
    "r": r_info,
    "stage2": {
        "isolated_library": rlib,
        "description_sha256": stage2_description_sha,
        "source_subdirectory": "stage2",
    },
    "activation_file": activate,
}
with open(manifest_path, "w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

rm -f "$OUT/python-info.json" "$OUT/r-info.json"
printf 'bootstrap: complete; manifest=%s activation=%s\n' "$OUT/environment.json" "$ACTIVATE"
