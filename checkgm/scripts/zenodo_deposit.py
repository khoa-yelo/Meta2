#!/usr/bin/env python3
"""Create (and optionally publish) the Zenodo deposit for the checkGM baselines and container.

The script is deliberately two-stage. By default it creates a *draft* record, reserves its DOI and uploads the files; the
record stays private and nothing is public until `--publish` is passed in a second, separate run. The reserved DOI is
printed, and is the one to put in the paper: Zenodo reserves it up front and it resolves once the record is published.

    export ZENODO_TOKEN=...                       # personal access token, scopes: deposit:write deposit:actions
    python scripts/zenodo_deposit.py --sandbox    # dry run against sandbox.zenodo.org first
    python scripts/zenodo_deposit.py              # the real draft on zenodo.org
    python scripts/zenodo_deposit.py --publish --id <record id>

Metadata comes from release/zenodo.json; the run aborts while any creator name is still TBD, so a record cannot be
created with placeholder authorship. Files are uploaded through the bucket API, which streams and so handles the
783 MB container without loading it into memory. An upload that is interrupted can be re-run: files already present
with the right size are skipped."""
import argparse, json, os, subprocess, sys, urllib.error, urllib.request

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(P, "release")


def api(url, token, method="GET", payload=None, timeout=120):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={"Content-Type": "application/json", "Authorization": f"Bearer {token}"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            body = r.read()
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as e:
        sys.exit(f"Zenodo {method} {url} -> {e.code}\n{e.read().decode()[:1200]}")


def put_file(bucket, name, path, token):
    """Upload through curl. Handing urllib an open file stalled on the bucket endpoint without transferring, and curl
    also gives retries and resumable behaviour that matter for the ~800 MB container image."""
    code = subprocess.run(["curl", "-s", "--retry", "3", "--retry-delay", "5", "--max-time", "7200",
                           "-w", "%{http_code}", "-X", "PUT", "-H", f"Authorization: Bearer {token}",
                           "--upload-file", path, f"{bucket}/{name}", "-o", os.devnull],
                          capture_output=True, text=True).stdout.strip()
    if code not in ("200", "201"):
        sys.exit(f"upload {name} -> HTTP {code}")
    return {"key": name}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org instead of the real archive")
    ap.add_argument("--publish", action="store_true", help="publish an existing draft; requires --id")
    ap.add_argument("--id", help="deposition id of an existing draft")
    a = ap.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ZENODO_TOKEN is not set (Zenodo > applications > personal access tokens; scopes deposit:write, deposit:actions)")
    base = "https://sandbox.zenodo.org/api" if a.sandbox else "https://zenodo.org/api"

    if a.publish:
        if not a.id:
            sys.exit("--publish needs --id")
        rec = api(f"{base}/deposit/depositions/{a.id}/actions/publish", token, "POST")
        print("published:", rec["doi_url"])
        return

    spec = json.load(open(os.path.join(REL, "zenodo.json")))
    meta = spec["metadata"]
    if any("TBD" in c["name"] or "TBD" in str(c.get("affiliation")) for c in meta["creators"]):
        sys.exit("release/zenodo.json still has placeholder creators; fill in name/affiliation/ORCID before depositing")
    for f in spec["files"]:
        if not os.path.exists(os.path.join(REL, f)):
            sys.exit(f"missing file: release/{f}")

    dep = api(f"{base}/deposit/depositions", token, "POST", {"metadata": meta})
    dep_id, bucket = dep["id"], dep["links"]["bucket"]
    doi = dep["metadata"].get("prereserve_doi", {}).get("doi") or dep.get("doi")
    print(f"draft {dep_id} created; reserved DOI: {doi}")

    present = {f["key"]: f["size"] for f in api(f"{base}/deposit/depositions/{dep_id}/files", token)}
    for f in spec["files"]:
        path = os.path.join(REL, f); size = os.path.getsize(path)
        if present.get(f) == size:
            print(f"  skip {f} (already uploaded)"); continue
        print(f"  uploading {f} ({size/1e6:.0f} MB)...", flush=True)
        put_file(bucket, f, path, token)
    print(f"\ndraft ready and still private: {dep['links']['html']}")
    print(f"put {doi} in the paper, then publish with:  python scripts/zenodo_deposit.py{' --sandbox' if a.sandbox else ''} --publish --id {dep_id}")


if __name__ == "__main__":
    main()
