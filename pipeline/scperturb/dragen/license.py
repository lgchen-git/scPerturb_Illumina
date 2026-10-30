#!/usr/bin/env python3
"""Check Illumina license service and credential without exposing the token."""

from __future__ import annotations

import argparse
import base64
import json
import stat
import sys
import urllib.error
import urllib.request
from pathlib import Path

from scperturb.common import load_json


def request_json(url: str, authorization: str | None = None) -> tuple[int, dict]:
    request = urllib.request.Request(
        url, headers={"Accept": "application/json", "Content-Type": "application/json"}
    )
    if authorization:
        request.add_header("Authorization", authorization)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with opener.open(request, timeout=30) as response:
            status, body = response.status, response.read()
    except urllib.error.HTTPError as exc:
        status, body = exc.code, exc.read()
    try:
        payload = json.loads(body)
    except json.JSONDecodeError:
        payload = {"message": "non-JSON response"}
    return status, payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--license-check", choices=("strict", "product-validation"), required=True)
    args = parser.parse_args()
    config = load_json(args.config)
    credential = Path(config["dragen"]["credential"])
    if not credential.is_file() or stat.S_IMODE(credential.stat().st_mode) != 0o600:
        raise PermissionError("DRAGEN credential must exist with mode 600")
    values = {}
    for line in credential.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    username, token = values.get("credentials-1"), values.get("credentials-2")
    if username != "IlluminaPlatform" or not token:
        raise ValueError("credential format is invalid; expected credentials-1 and credentials-2")
    health_status, health = request_json("https://license.dragen.illumina.com/healthcheck/version")
    print(f"healthcheck_status={health_status}")
    if "version" in health:
        print(f"license_server_version={health['version']}")
    if health_status != 200:
        raise RuntimeError("license server healthcheck failed")
    basic = base64.b64encode(f"{username}:{token}".encode()).decode()
    quota_status, quota = request_json("https://license.dragen.illumina.com/api/v2/query_quota", "Basic " + basic)
    print(f"quota_status={quota_status}")
    print(f"quota_message={quota.get('message', '')}")
    licenses = quota.get("licenses") if isinstance(quota.get("licenses"), list) else []
    print(f"license_count={len(licenses)}")
    if quota_status == 200 and licenses:
        print("LICENSE_PREFLIGHT=PASS")
        return 0
    if args.license_check == "product-validation" and quota_status == 403 and "Authentication Failed" in str(quota.get("message", "")):
        print("LICENSE_PREFLIGHT=INCONCLUSIVE_QUOTA_403")
        return 0
    raise RuntimeError("license quota authentication failed or no active license was returned")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
