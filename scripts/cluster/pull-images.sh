#!/usr/bin/env bash
# Run on a login/transfer node with network access, before sbatch.
set -euo pipefail
if [[ $# != 2 || ! $1 =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-][A-Za-z0-9.-]+)?$ ]]; then
    echo "Usage: bash scripts/cluster/pull-images.sh vX.Y.Z NEW_IMAGE_DIRECTORY" >&2
    exit 2
fi
runtime=${SEGMETER_RUNTIME:-singularity}
case "$runtime" in singularity|apptainer) ;; *) echo "Unsupported runtime: $runtime" >&2; exit 2 ;; esac
command -v "$runtime" >/dev/null
mkdir "$2"
image_dir=$(cd "$2" && pwd)
for container in others giggle rust-tools; do
    "$runtime" pull "$image_dir/$container.sif" "docker://yanglabinfo/segmeter:$container-$1"
done
echo "Images ready in $image_dir"
