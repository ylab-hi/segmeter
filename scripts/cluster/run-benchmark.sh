#!/usr/bin/env bash
# Prepare on a login node and submit the offline benchmark to Slurm.
set -euo pipefail
usage() {
    cat <<'EOF'
Usage: bash scripts/cluster/run-benchmark.sh --image-tag vX.Y.Z --work-dir /shared/path [options]
  --module NAME          Singularity module name (default: singularity; none skips loading)
  --runtime NAME         singularity or apptainer (default: singularity)
  --dataset NAME         simulated or zenodo (default: simulated)
  --account NAME         Slurm account
  --partition NAME       Slurm partition
  --sizes LIST           Optional size override, e.g. 1K for a setup check
  --mem SIZE             Slurm memory (default: 32G)
  --time LIMIT           Slurm time limit (default: 1-00:00:00)
  --exclusive            Request an exclusive node
Images and Zenodo data are cached under WORK_DIR. Each submission gets a new results directory.
EOF
}
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image_tag='' work_dir='' sizes='' dataset=simulated
export SEGMETER_MODULE=singularity SEGMETER_RUNTIME=singularity
slurm_options=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --exclusive) slurm_options+=(--exclusive); shift; continue ;;
        --image-tag|--work-dir|--module|--runtime|--dataset|--account|--partition|--sizes|--mem|--time)
            [[ $# -ge 2 ]] || { usage >&2; exit 2; } ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    case "$1" in
        --image-tag) image_tag=$2 ;;
        --work-dir) work_dir=$2 ;;
        --module) SEGMETER_MODULE=$2 ;;
        --runtime) SEGMETER_RUNTIME=$2 ;;
        --dataset) dataset=$2 ;;
        --sizes) sizes=$2 ;;
        --account|--partition|--mem|--time) slurm_options+=("$1=$2") ;;
    esac
    shift 2
done
[[ -n $work_dir && $image_tag =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-][A-Za-z0-9.-]+)?$ ]] || { usage >&2; exit 2; }
case "$dataset" in simulated|zenodo) ;; *) echo "Invalid dataset: $dataset" >&2; exit 2 ;; esac
case "$SEGMETER_RUNTIME" in singularity|apptainer) ;; *) echo "Invalid runtime" >&2; exit 2 ;; esac
# shellcheck source=scripts/cluster/load-runtime.sh
source "$repo/scripts/cluster/load-runtime.sh"
load_segmeter_runtime
command -v python3 >/dev/null
command -v sbatch >/dev/null
command -v flock >/dev/null
mkdir -p "$work_dir"
work_dir=$(cd "$work_dir" && pwd)
image_dir="$work_dir/images/$image_tag"
mkdir -p "$image_dir" "$work_dir/results"
# A lock prevents simultaneous launchers from pulling/extracting into the same cache.
exec 9>"$work_dir/.prepare.lock"
flock -n 9 || { echo "Another launcher is preparing $work_dir" >&2; exit 1; }
for container in others giggle rust-tools; do
    if [[ ! -f $image_dir/$container.sif ]]; then
        "$SEGMETER_RUNTIME" pull --force "$image_dir/$container.sif.part" "docker://yanglabinfo/segmeter:$container-$image_tag"
        mv "$image_dir/$container.sif.part" "$image_dir/$container.sif"
    fi
done
benchmark_options=()
if [[ $dataset == zenodo ]]; then
    python3 "$repo/scripts/cluster/prepare_zenodo.py" --directory "$work_dir/zenodo-14880992"
    benchmark_options+=(--dataset-dir "$work_dir/zenodo-14880992/simdata")
fi
[[ -z $sizes ]] || benchmark_options+=(--sizes "$sizes")
flock -u 9
stamp=$(date -u +%Y%m%dT%H%M%SZ)
run_dir=$(mktemp -d "$work_dir/results/$dataset-$stamp-XXXXXX")
# Freeze the checkout now, so later edits while queued cannot change this run.
mkdir "$run_dir/checkout"
cp -R "$repo/segmeter" "$repo/containers" "$repo/scripts" "$run_dir/checkout/"
git -C "$repo" rev-parse HEAD > "$run_dir/checkout/source-commit.txt"
cd "$run_dir/checkout"
job_id=$(sbatch --parsable "${slurm_options[@]}" \
    --output="$run_dir/slurm-%j.log" \
    scripts/cluster/full-benchmark.sbatch "$image_dir" "$run_dir/benchmark" "${benchmark_options[@]}")
printf 'Submitted Slurm job %s\nResults: %s/benchmark\nLog: %s/slurm-%s.log\n' "$job_id" "$run_dir" "$run_dir" "${job_id%%;*}"
