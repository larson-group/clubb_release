#!/usr/bin/env bash
# Set up the dedicated Python environment and run CLUBB JAX.
#
# This is Bash so it can be invoked directly and install a supported Python when
# one is not available. It owns all JAX environment setup, then replaces itself
# with the CLUBB JAX Python process. Callers do not need to install dependencies,
# activate a virtualenv, or know how the environment is managed.
set -euo pipefail

# Repository paths and overridable locations. Keeping managed tools separate from
# the virtualenv lets us replace a stale environment without losing uv or Python.
readonly uv_version="0.11.32"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/.." && pwd)"
accelerator="${CLUBB_JAX_ACCELERATOR:-cpu}"
accelerator="${accelerator,,}"
case "$accelerator" in
  cpu)
    requirements="$script_dir/requirements.txt"
    default_venv="$repo_root/.venv-jax"
    jax_platform="cpu"
    ;;
  cuda13)
    requirements="$script_dir/requirements-cuda13.txt"
    default_venv="$repo_root/.venv-jax-cuda13"
    # JAX_PLATFORMS expects the concrete plugin name. "gpu" also probes ROCm.
    jax_platform="cuda"
    ;;
  *)
    echo "ERROR: CLUBB_JAX_ACCELERATOR must be 'cpu' or 'cuda13'; got: $accelerator" >&2
    exit 1
    ;;
esac
export CLUBB_JAX_ACCELERATOR="$accelerator"
export JAX_PLATFORMS="$jax_platform"
venv_dir="${CLUBB_JAX_VENV:-$default_venv}"
tools_dir="${CLUBB_JAX_TOOLS_DIR:-$repo_root/.clubb-jax-tools}"
[[ "$venv_dir" == /* ]] || venv_dir="$repo_root/$venv_dir"
[[ "$tools_dir" == /* ]] || tools_dir="$repo_root/$tools_dir"

# Command-line help and small shared utilities.
usage() {
  cat <<'EOF'
Usage: clubb_jax/run_jax_wrapper.sh [namelist-path] [driver-args...]
       clubb_jax/run_jax_wrapper.sh --init_env

Creates or updates the JAX environment, then runs the CLUBB JAX standalone.
Missing copies of uv and Python are downloaded automatically.

Environment:
  CLUBB_JAX_ACCELERATOR  Runtime backend: cpu (default) or cuda13
  CLUBB_JAX_VENV      Virtualenv path (default: .venv-jax or .venv-jax-cuda13)
  CLUBB_JAX_TOOLS_DIR Managed uv/Python path (default: .clubb-jax-tools)
  PYTHON               Python used when creating a new virtualenv
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Python/JAX compatibility policy. Python 3.12+ gets JAX 0.11.0; Python 3.11
# remains a supported fallback with JAX 0.10.0.
python_is_supported() {
  "$1" -c 'import sys; raise SystemExit(sys.version_info < (3, 11))' >/dev/null 2>&1
}

python_version() {
  "$1" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")'
}

jax_version_for() {
  if "$1" -c 'import sys; raise SystemExit(sys.version_info < (3, 12))'; then
    echo 0.11.0
  else
    echo 0.10.0
  fi
}

# Honor an explicit PYTHON first. Otherwise prefer the newest supported Python
# already on PATH so machines with a suitable interpreter avoid a download.
find_python() {
  if [[ -n "${PYTHON:-}" ]]; then
    command -v "$PYTHON" >/dev/null 2>&1 || die "PYTHON does not exist: $PYTHON"
    python_is_supported "$PYTHON" ||
      die "JAX 0.10.0 or newer requires Python 3.11+; selected: $($PYTHON --version 2>&1)"
    echo "$PYTHON"
    return
  fi

  local candidate
  for candidate in python3.14 python3.13 python3.12 python3.11 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 && python_is_supported "$candidate"; then
      command -v "$candidate"
      return
    fi
  done
}

# Reuse uv when available, or install a pinned copy without changing PATH or shell
# profiles. Pinning makes first-run behavior stable instead of tracking latest uv.
ensure_uv() {
  local local_uv="$tools_dir/bin/uv"
  if [[ -x "$local_uv" ]]; then
    uv_cmd="$local_uv"
    return
  fi
  if command -v uv >/dev/null 2>&1; then
    uv_cmd="$(command -v uv)"
    return
  fi

  local installer="https://astral.sh/uv/$uv_version/install.sh"
  mkdir -p "$tools_dir/bin"
  echo "==> Installing uv $uv_version in $tools_dir/bin"
  if command -v curl >/dev/null 2>&1; then
    curl --proto '=https' --tlsv1.2 -LsSf "$installer" |
      env UV_UNMANAGED_INSTALL="$tools_dir/bin" sh
  elif command -v wget >/dev/null 2>&1; then
    wget -qO- "$installer" |
      env UV_UNMANAGED_INSTALL="$tools_dir/bin" sh
  else
    die "Automatic setup requires curl or wget to download uv."
  fi

  [[ -x "$local_uv" ]] || die "uv installation did not create $local_uv"
  uv_cmd="$local_uv"
  echo "==> uv is ready: $($uv_cmd --version)"
}

# Parse launcher-only options before forwarding all remaining arguments unchanged
# to the standalone driver.
case "${1:-}" in
  --launcher-help) usage; exit 0 ;;
  --init_env) init_env_only=true; shift ;;
  *) init_env_only=false ;;
esac

# Keep uv's downloaded Python and cache inside the ignored tools directory. This
# avoids administrator access and leaves the system Python untouched.
[[ -f "$requirements" ]] || die "Missing requirements file: $requirements"
mkdir -p "$tools_dir"
export UV_CACHE_DIR="$tools_dir/cache"
export UV_PYTHON_INSTALL_DIR="$tools_dir/python"

# Parallel comparison runs share one environment. Locking outside the virtualenv
# prevents two processes from creating, clearing, or installing into it at once.
exec 9>"$tools_dir/environment.lock"
flock 9

# Reuse a healthy environment. If the default environment is missing or stale,
# create it from an installed Python or let uv download Python 3.12. Custom venvs
# are never cleared automatically because they may contain user-managed state.
venv_python="$venv_dir/bin/python"
if [[ -x "$venv_python" ]] && python_is_supported "$venv_python"; then
  echo "==> Reusing JAX virtualenv: $venv_dir (Python $(python_version "$venv_python"))"
else
  if [[ -e "$venv_dir" && -n "${CLUBB_JAX_VENV:-}" ]]; then
    die "Custom virtualenv is incomplete or uses Python older than 3.11: $venv_dir"
  fi

  python_cmd="$(find_python)"
  if [[ -n "$python_cmd" ]]; then
    python_spec="$python_cmd"
    echo "==> Creating JAX virtualenv with existing Python $(python_version "$python_cmd")"
  else
    python_spec="3.12"
    echo "==> No Python 3.11+ found; uv will download Python 3.12"
  fi

  ensure_uv
  if [[ -e "$venv_dir" ]]; then
    echo "==> Replacing unsupported JAX virtualenv: $venv_dir"
    "$uv_cmd" venv --clear --python "$python_spec" "$venv_dir"
  else
    "$uv_cmd" venv --python "$python_spec" "$venv_dir"
  fi
fi

# A content hash avoids running the installer on every launch. Package versions
# are also checked so a stale or manually altered environment repairs itself.
required_jax="$(jax_version_for "$venv_python")"
requirements_hash="$($venv_python -c 'import hashlib, pathlib, sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$requirements")"
requirements_stamp="$venv_dir/.clubb-jax-requirements.sha256"
installed_hash=""
[[ ! -f "$requirements_stamp" ]] || installed_hash="$(<"$requirements_stamp")"

packages_are_ready() {
  "$venv_python" -c '
from importlib.metadata import version
import sys

assert version("jax") == sys.argv[1]
assert version("jaxlib") == sys.argv[1]
for package in ("netCDF4", "pytest", "tabulate"):
    version(package)
if sys.argv[2] == "cuda13":
    assert version("jax-cuda13-plugin") == sys.argv[1]
    assert version("jax-cuda13-pjrt") == sys.argv[1]
' "$required_jax" "$accelerator" >/dev/null 2>&1
}

environment_is_ready() {
  [[ "$requirements_hash" == "$installed_hash" ]] && packages_are_ready
}

if ! environment_is_ready; then
  ensure_uv
  echo "==> Installing JAX $required_jax ($accelerator) and requirements from $requirements"
  "$uv_cmd" pip install --python "$venv_python" -r "$requirements"
  packages_are_ready || die "JAX $accelerator packages failed validation after installation."
  printf '%s\n' "$requirements_hash" >"$requirements_stamp"
fi

# Installation is complete, so other waiting runs may safely inspect and reuse it.
flock -u 9

# --init_env stops after reporting the resolved runtime. Normal invocations use
# exec so signals and the driver's exit status pass directly back to run_scm.py.
if [[ "$init_env_only" == true ]]; then
  "$venv_python" -c 'import jax, jaxlib; print(f"JAX environment ready: jax={jax.__version__} jaxlib={jaxlib.__version__} backend={jax.default_backend()} devices={jax.devices()}")'
  exit
fi

cd "$repo_root"
exec "$venv_python" -m clubb_jax.src.clubb_standalone "$@"
