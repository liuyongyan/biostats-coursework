"""
Query GDC REST API for TCGA DNA methylation sample counts per project.

For each TCGA project, we want:
- N unique cases (patients) with 450K methylation data
- N unique cases with 27K methylation data
- N files on each platform (per sample-type: primary tumor vs solid tissue normal)
- Paired-patient count (patients with BOTH primary tumor AND adjacent normal 450K)
- Estimated total file size for 450K methylation beta-value files

Output: prints a markdown table and writes JSON to gdc_methylation_counts.json.
"""

import json
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

BASE = "https://api.gdc.cancer.gov"


def post_json(endpoint, payload):
    req = urllib.request.Request(
        f"{BASE}/{endpoint}",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def get_json(endpoint, params):
    q = urllib.parse.urlencode(params)
    with urllib.request.urlopen(f"{BASE}/{endpoint}?{q}", timeout=120) as r:
        return json.loads(r.read())


def fetch_files(platform):
    """
    Fetch all TCGA methylation beta-value file records for a given platform.
    Paginate.
    """
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.program.name", "value": ["TCGA"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["DNA Methylation"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
            {"op": "in", "content": {"field": "platform", "value": [platform]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    fields = ",".join([
        "file_id",
        "file_size",
        "cases.project.project_id",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
    ])
    all_hits = []
    page_size = 1000
    frm = 0
    while True:
        payload = {
            "filters": filters,
            "fields": fields,
            "format": "JSON",
            "size": str(page_size),
            "from": str(frm),
        }
        resp = post_json("files", payload)
        hits = resp["data"]["hits"]
        all_hits.extend(hits)
        total = resp["data"]["pagination"]["total"]
        frm += page_size
        if frm >= total:
            break
    return all_hits


def summarize(hits_450, hits_27):
    # project -> { cases_450: set, cases_27: set, tumor_cases, normal_cases, size_bytes, files_n }
    proj = defaultdict(lambda: {
        "cases_450": set(),
        "cases_27": set(),
        "tumor_cases_450": set(),      # case IDs with primary tumor on 450K
        "normal_cases_450": set(),      # case IDs with solid tissue normal on 450K
        "metastatic_cases_450": set(),
        "files_n_450": 0,
        "files_n_27": 0,
        "size_bytes_450": 0,
        "size_bytes_27": 0,
    })

    def ingest(hits, platform_tag):
        for h in hits:
            case = h["cases"][0]
            pid = case["project"]["project_id"]
            sub = case["submitter_id"]
            size = h.get("file_size", 0) or 0
            samples = case.get("samples", [])
            sample_types = [s.get("sample_type", "") for s in samples]
            p = proj[pid]
            if platform_tag == "450":
                p["cases_450"].add(sub)
                p["files_n_450"] += 1
                p["size_bytes_450"] += size
                # Sample-type classification at file level
                for st in sample_types:
                    if "Primary Tumor" in st or "Primary Blood Derived" in st:
                        p["tumor_cases_450"].add(sub)
                    elif "Solid Tissue Normal" in st or "Blood Derived Normal" in st:
                        p["normal_cases_450"].add(sub)
                    elif "Metastatic" in st:
                        p["metastatic_cases_450"].add(sub)
            else:
                p["cases_27"].add(sub)
                p["files_n_27"] += 1
                p["size_bytes_27"] += size

    ingest(hits_450, "450")
    ingest(hits_27, "27")
    return proj


def main():
    print("Querying GDC for 450K methylation files ...", file=sys.stderr)
    hits_450 = fetch_files("Illumina Human Methylation 450")
    print(f"  got {len(hits_450)} files", file=sys.stderr)

    print("Querying GDC for 27K methylation files ...", file=sys.stderr)
    hits_27 = fetch_files("Illumina Human Methylation 27")
    print(f"  got {len(hits_27)} files", file=sys.stderr)

    proj = summarize(hits_450, hits_27)

    rows = []
    for pid, d in proj.items():
        paired = d["tumor_cases_450"] & d["normal_cases_450"]
        rows.append({
            "project": pid,
            "cases_450": len(d["cases_450"]),
            "tumor_450": len(d["tumor_cases_450"]),
            "normal_450": len(d["normal_cases_450"]),
            "metastatic_450": len(d["metastatic_cases_450"]),
            "paired_tumor_normal_450": len(paired),
            "files_450": d["files_n_450"],
            "size_gb_450": round(d["size_bytes_450"] / 1e9, 2),
            "cases_27": len(d["cases_27"]),
            "size_gb_27": round(d["size_bytes_27"] / 1e9, 2),
        })

    # Sort by number of paired 450K samples descending, then by tumor count
    rows.sort(key=lambda r: (-r["paired_tumor_normal_450"], -r["tumor_450"]))

    # Print markdown table
    print("| Project | Cases 450K | Tumor 450K | Normal 450K | Metastatic 450K | **Paired 450K** | Files 450K | Size 450K (GB) | Cases 27K | Size 27K (GB) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['project']} | {r['cases_450']} | {r['tumor_450']} | {r['normal_450']} | {r['metastatic_450']} | **{r['paired_tumor_normal_450']}** | {r['files_450']} | {r['size_gb_450']} | {r['cases_27']} | {r['size_gb_27']} |")

    # Totals
    total_files = sum(r["files_450"] for r in rows)
    total_gb = sum(r["size_gb_450"] for r in rows)
    print()
    print(f"Total 450K files: {total_files}  |  Total 450K size: {total_gb:.1f} GB")

    with open("gdc_methylation_counts.json", "w") as f:
        json.dump(rows, f, indent=2)


if __name__ == "__main__":
    main()
