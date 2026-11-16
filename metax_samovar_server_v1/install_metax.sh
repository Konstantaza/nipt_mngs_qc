#!/usr/bin/env bash
set -euo pipefail
# Run from the nipt_mngs_qc repository, with its main environment active.
test -f nextflow.config
test -f modules/local/annotate_metax.nf
repo="$(pwd -P)"
prefix="$(dirname "$repo")/metax_env_2_0_3"
kit="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
conda_bin="${CONDA_EXE:-$(command -v conda)}"
if [[ -e "$prefix" ]]; then
    test -x "$prefix/bin/metax"
    [[ "$("$prefix/bin/metax" --version)" == "metax 2.0.3" ]]
    printf 'Reuse verified MetaX 2.0.3: %s\n' "$prefix"
else
    "$conda_bin" create --yes --prefix "$prefix" \
        --override-channels -c zldeng -c conda-forge \
        metax=2.0.3 python=3.14.7 2>&1 | tee "$kit/install.log"
fi
"$prefix/bin/metax" --version
"$conda_bin" list --prefix "$prefix" --explicit > "$kit/metax_environment.explicit.txt"
printf 'Main Python is unchanged: %s\n' "$(command -v python)"
