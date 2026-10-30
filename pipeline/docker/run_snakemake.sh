#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
config_rel=config/run.yaml
cores=${SC_PERTURB_CORES:-8}
mode=full
dry_run=false
while (($#)); do
  case "$1" in
    --config) config_rel=$2; shift 2 ;;
    --cores) cores=$2; shift 2 ;;
    --pilot-only) [[ "$mode" == full ]] || { echo "run modes are mutually exclusive" >&2; exit 2; }; mode=pilot; shift ;;
    --report-only) [[ "$mode" == full ]] || { echo "run modes are mutually exclusive" >&2; exit 2; }; mode=report; shift ;;
    --review-only) [[ "$mode" == full ]] || { echo "run modes are mutually exclusive" >&2; exit 2; }; mode=review; shift ;;
    --through-dragen) [[ "$mode" == full ]] || { echo "run modes are mutually exclusive" >&2; exit 2; }; mode=dragen; shift ;;
    --prepare-icm) [[ "$mode" == full ]] || { echo "run modes are mutually exclusive" >&2; exit 2; }; mode=icm; shift ;;
    -h|--help) echo "Usage: $0 [--config config/run.yaml] [--cores 8] [--through-dragen | --prepare-icm | --pilot-only | --report-only | --review-only] [--dry-run]"; exit 0 ;;
    --dry-run) dry_run=true; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
config_host="$project_root/$config_rel"
export PYTHONPATH="$project_root/pipeline${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
image_tag=${SC_PERTURB_IMAGE:-scperturb-crispr:dragen4.5.4-1-v1.0.0}
[[ -f "$config_host" ]] || { echo "missing run config: $config_host" >&2; exit 2; }
[[ "$cores" =~ ^[1-9][0-9]*$ ]] || { echo "--cores must be a positive integer" >&2; exit 2; }
if [[ "$mode" == icm ]]; then
  icm_args=(prepare-icm --project "$project_root" --config "$config_host")
  if [[ "$dry_run" == true ]]; then icm_args+=(--dry-run); fi
  python3 -m scperturb.runner "${icm_args[@]}"
  exit 0
fi
if [[ "$mode" == report || "$mode" == review ]]; then
  python3 -m scperturb.runner reports-ready --project "$project_root" --config "$config_host"
  if [[ "$dry_run" == true ]]; then
    echo "RUN_DRY_RUN mode=$mode upstream_stages=none network=none config=$config_host"
    exit 0
  fi
else
  reference_root=${SC_PERTURB_HG38_ROOT:-}
  if [[ -z "$reference_root" && -f "$project_root/config/reference_path.txt" ]]; then
    IFS= read -r reference_root < "$project_root/config/reference_path.txt" || true
  fi
  [[ -n "$reference_root" ]] || { echo "set SC_PERTURB_HG38_ROOT or add validated hg38 root to config/reference_path.txt" >&2; exit 2; }
  credential=${DRAGEN_CREDENTIAL:-}
  if [[ -z "$credential" && -f "$project_root/config/credential_path.txt" ]]; then
    IFS= read -r credential < "$project_root/config/credential_path.txt" || true
  fi
  [[ -n "$credential" ]] || { echo "set DRAGEN_CREDENTIAL or add the mode-600 license path to config/credential_path.txt" >&2; exit 2; }
  python3 -m scperturb.runner validate --project "$project_root" --config "$config_host" --reference "$reference_root" --credential "$credential"
fi
docker image inspect "$image_tag" >/dev/null
run_id=$(date -u +%Y%m%dT%H%M%SZ)
git_head=$(git -C "$project_root" rev-parse --verify HEAD 2>/dev/null || printf 'NO_COMMIT')
git_status=$(git -C "$project_root" status --porcelain=v1 --untracked-files=normal 2>/dev/null | base64 -w0)
mkdir -p "$project_root/log" "$project_root/outputs" "$project_root/.snakemake" "$project_root/tmp/runner/$run_id/home"
log_file="$project_root/log/snakemake_$run_id.log"
container_config="/project/$config_rel"
if [[ "$mode" == pilot ]]; then
  pilot_config="$project_root/tmp/runner/$run_id/pilot_config.json"
  python3 - "$config_host" "$pilot_config" "$run_id" <<'PY'
from pathlib import Path
import json, sys
source, destination, run_id = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
config = json.loads(source.read_text(encoding="utf-8"))
root = f"/project/tmp/native_license_pilot_{run_id}"
config["results_root"] = f"{root}/outputs"
config["log_root"] = f"{root}/log"
config["tmp_root"] = f"{root}/tmp"
config["pilot_only"] = True
destination.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
  container_config="/project/tmp/runner/$run_id/pilot_config.json"
fi
args=(--snakefile pipeline/Snakefile --configfile "$container_config" --cores "$cores" --resources "mem_mb=${SC_PERTURB_MEMORY_MB:-32768}" --rerun-incomplete --printshellcmds --show-failed-logs)
if [[ "$mode" == dragen ]]; then args+=(--until dragen_scrna_crispr); fi
if [[ "$dry_run" == true ]]; then args+=(--dry-run); fi
network=host
entrypoint_args=()
mount_args=()
if [[ "$mode" == report || "$mode" == review ]]; then
  network=none
  entrypoint_args=(--entrypoint /opt/conda/envs/illumina_scrna/bin/python)
  subcommand=refresh
  if [[ "$mode" == review ]]; then subcommand=review; fi
  args=(-m scperturb.runner "$subcommand" --project /project --config "$container_config")
else
  mount_args+=(--mount "type=bind,source=$reference_root,target=/reference/hg38,readonly"
               --mount "type=bind,source=$credential,target=/run/secrets/dragen-credential,readonly")
fi
docker_args=(run --rm --network "$network" --ulimit nofile=65535:65535 --user "$(id -u):$(id -g)" --workdir /project
  --env "HOME=/project/tmp/runner/$run_id/home" --env "TMPDIR=/project/tmp/runner/$run_id" --env "TMP=/project/tmp/runner/$run_id" --env "TEMP=/project/tmp/runner/$run_id"
  --env "PYTHONPATH=/project/pipeline" --env "PYTHONDONTWRITEBYTECODE=1"
  --env "NUMBA_CACHE_DIR=/project/tmp/runner/$run_id/numba"
  --env "XDG_CACHE_HOME=/project/tmp/runner/$run_id/cache" --env "MPLCONFIGDIR=/project/tmp/runner/$run_id/mpl"
  --mount "type=bind,source=$project_root/config,target=/project/config,readonly"
  --mount "type=bind,source=$project_root/inputs,target=/project/inputs,readonly"
  --mount "type=bind,source=$project_root/reference,target=/project/reference,readonly"
  --mount "type=bind,source=$project_root/pipeline,target=/project/pipeline,readonly"
  --mount "type=bind,source=$project_root/.snakemake,target=/project/.snakemake"
  --mount "type=bind,source=$project_root/log,target=/project/log"
  --mount "type=bind,source=$project_root/outputs,target=/project/outputs"
  --mount "type=bind,source=$project_root/tmp,target=/project/tmp"
  --env "SC_PERTURB_GIT_HEAD=$git_head" --env "SC_PERTURB_GIT_STATUS_BASE64=$git_status"
  "${mount_args[@]}" "${entrypoint_args[@]}" "$image_tag" "${args[@]}")
set +e
docker "${docker_args[@]}" >"$log_file" 2>&1
status=$?
set -e
if ((status != 0)); then
  echo "Workflow failed (mode=$mode exit=$status); log: $log_file" >&2
  tail -n 50 "$log_file" >&2
  exit "$status"
fi
if [[ "$dry_run" == true ]]; then
  echo "RUN_DRY_RUN log=$log_file"
elif [[ "$mode" == pilot ]]; then
  python3 - "$project_root" "$run_id" <<'PY'
from pathlib import Path
import json, sys
root=Path(sys.argv[1]); run_id=sys.argv[2]
base=root/"tmp"/f"native_license_pilot_{run_id}"/"outputs"
markers=list(base.glob("*/02_dragen_scrna_crispr/05_metadata/stage_complete.json"))
if not markers: raise SystemExit("Snakemake returned success without a native DRAGEN pilot marker")
for path in markers:
    doc=json.loads(path.read_text(encoding="utf-8"))
    if doc.get("status")!="complete" or doc.get("interface")!="02-dragen-license-pilot-v1": raise SystemExit(f"invalid pilot marker: {path}")
print(f"PILOT_COMPLETE samples={len(markers)} scope=100000_reads")
PY
elif [[ "$mode" == dragen ]]; then
  python3 -m scperturb.runner prepare-icm --project "$project_root" --config "$config_host"
else
  python3 -m scperturb.runner verify --project "$project_root" --config "$config_host"
fi
if [[ "$dry_run" != true ]]; then echo "RUN_COMPLETE mode=$mode log=$log_file"; fi
