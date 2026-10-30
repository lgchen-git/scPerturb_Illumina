#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image_tag=${SC_PERTURB_IMAGE:-scperturb-crispr:dragen4.5.4-1-v1.0.0}
archive_path=${1:-}
read -r base_image expected_id < "$project_root/pipeline/base-image.lock"
actual_id=$(docker image inspect --format '{{.Id}}' "$base_image")
if [[ "$actual_id" != "$expected_id" ]]; then
  echo "base image mismatch: $base_image has $actual_id; expected $expected_id" >&2
  exit 2
fi
docker build --pull=false --network=host --file "$project_root/pipeline/Dockerfile" --tag "$image_tag" "$project_root/pipeline"
if [[ -n "$archive_path" ]]; then
  mkdir -p "$(dirname "$archive_path")"
  docker save --output "$archive_path" "$image_tag"
fi
