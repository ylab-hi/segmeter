#!/usr/bin/env bash
# Sourced by the launcher and batch job; site module name is configurable.
load_segmeter_runtime() {
    local init_file
    if [[ ${SEGMETER_MODULE:-singularity} != none || -n ${SEGMETER_PYTHON_MODULE:-} ]]; then
        if ! type module >/dev/null 2>&1; then
            for init_file in /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
                if [[ -r $init_file ]]; then
                    # shellcheck disable=SC1090
                    source "$init_file"
                    break
                fi
            done
        fi
        if ! type module >/dev/null 2>&1; then
            echo "Module command unavailable. Initialize your site's modules, or use --module none if Singularity is already on PATH." >&2
            return 1
        fi
        if [[ ${SEGMETER_MODULE:-singularity} != none ]]; then
            module load "${SEGMETER_MODULE:-singularity}"
        fi
        if [[ -n ${SEGMETER_PYTHON_MODULE:-} ]]; then
            module load "$SEGMETER_PYTHON_MODULE"
        fi
    fi
    command -v "${SEGMETER_RUNTIME:-singularity}" >/dev/null
    command -v python3 >/dev/null
    python3 -c 'import sys; sys.exit("Host Python 3.8 or newer is required; use --python-module with your cluster Python module" if sys.version_info < (3, 8) else 0)'
}
