#!/usr/bin/env bash
# Separate generator Python; do not downgrade the pipeline's Python 3.14.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMOVAR="${1:-$HERE/../../../../samovar}"
SAMOVAR="$(cd "$SAMOVAR" && pwd)"
VENV="${2:-$SAMOVAR/.venv-nipt-toy}"
if [[ -e "$VENV" ]]; then
    echo "Refusing to overwrite existing environment: $VENV" >&2
    echo "Reuse its bin/python, or choose a new second argument." >&2
    exit 1
fi
command -v uv >/dev/null || {
    echo "Install uv first, e.g. conda install -c conda-forge uv" >&2
    exit 1
}
test -f "$SAMOVAR/src/samovar/iss_config.py"
uv venv --python 3.11 "$VENV"
uv pip install --python "$VENV/bin/python" -r "$HERE/samovar_requirements.lock.txt" -e "$SAMOVAR"
"$VENV/bin/python" -c 'import sys,importlib.metadata as m; assert sys.version_info[:2]==(3,11); print(sys.version); print({p:m.version(p) for p in ("samovar","insilicoseq","snakemake")})'
echo "Generator interpreter: $VENV/bin/python"
