#!/usr/bin/env python3
"""Refresh the bundled real ENA/NCBI lineage subset; normal runs are offline."""
import json
from pathlib import Path
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/ena/taxonomy/rest/"
HERE = Path(__file__).resolve().parent


def fetch_subset(destination=HERE/"taxonomy"):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    cache_file = destination/"ena_responses.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}

    def get(endpoint):
        url = BASE + endpoint
        if url not in cache:
            with urllib.request.urlopen(url, timeout=90) as r:
                cache[url] = json.load(r)
            cache_file.write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n")
        return cache[url]

    nodes = {"1": ("1", "no rank", "root")}
    for leaf in ("2886930", "10847", "4932", "562", "9606"):
        record = get("tax-id/" + leaf)
        lineage = [s.strip() for s in record["lineage"].split(";") if s.strip()]
        parent = "1"
        for depth, name in enumerate(lineage + [record["scientificName"]]):
            if name == record["scientificName"]:
                item = record
            else:
                candidates = get("scientific-name/" + urllib.parse.quote(name, safe=""))
                if isinstance(candidates, dict):
                    candidates = [candidates]
                exact = [x for x in candidates if x["scientificName"] == name and
                         [s.strip() for s in x.get("lineage", "").split(";") if s.strip()]
                         == lineage[:depth]]
                if len(exact) != 1:
                    raise ValueError(f"Ambiguous taxonomy name: {name}")
                item = exact[0]
            tid = str(item["taxId"])
            node = (parent, item["rank"], name)
            if tid in nodes and nodes[tid] != node:
                raise ValueError(f"Inconsistent lineage {tid}: {node} vs {nodes[tid]}")
            nodes[tid] = node
            parent = tid
    if nodes["2886930"][:2] != ("10847", "no rank") or nodes["10847"][1] != "species":
        raise ValueError("phiX174 taxonomy changed; review mapping before refreshing")
    with open(destination/"nodes.dmp", "w") as f:
        for tid, (parent, rank, name) in sorted(nodes.items(), key=lambda x: int(x[0])):
            f.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\n")
    with open(destination/"names.dmp", "w") as f:
        for tid, (parent, rank, name) in sorted(nodes.items(), key=lambda x: int(x[0])):
            f.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    print(f"Saved {len(nodes)} real taxonomy nodes to {destination}")


if __name__ == "__main__":
    fetch_subset()
