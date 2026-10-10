#!/usr/bin/env bash
# Site-specific account, partition and resource overrides belong in the sbatch command.
#SBATCH --job-name=segmeter-full
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=segmeter-full-%j.log

set -euo pipefail
if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch scripts/cluster/full-benchmark.sbatch IMAGE_DIRECTORY NEW_OUTPUT_DIRECTORY [benchmark arguments]" >&2
    exit 2
fi
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the segmeter checkout}"
runtime=${SEGMETER_RUNTIME:-singularity}
case "$runtime" in singularity|apptainer) ;; *) echo "Unsupported runtime: $runtime" >&2; exit 2 ;; esac
if [[ -n ${SEGMETER_MODULE:-} ]]; then
    # shellcheck source=scripts/cluster/load-runtime.sh
    source scripts/cluster/load-runtime.sh
    load_segmeter_runtime
fi
command -v "$runtime" >/dev/null
command -v python3 >/dev/null
image_dir=$1
output_dir=$2
shift 2
# All tools run sequentially in this single allocation. Images must be pre-pulled;
# the job requires neither network access, Docker, nor GitHub credentials.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/continuous_benchmark.py \
    --tier full --runtime "$runtime" --sif-dir "$image_dir" \
    --output "$output_dir" "$@"
echo "Benchmark complete. Publication files: $output_dir/publish/"
