#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${PIPELINE_PYTHON:-python3}" "$HERE/toy.py" "$@"
