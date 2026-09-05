"""Recheck Sarid's rational certificate with an independent implementation.

This program enumerates all 2^15 labelled six-vertex graphs itself. It does
not import Sarid's verifier, NetworkX, or a precomputed graph catalogue.
The mathematical interpretation of the flag squares is given in the paper.
"""

import argparse
import hashlib
import itertools
import json
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


def edge_pairs(n):
    return list(itertools.combinations(range(n), 2))


def adjacency(n, mask):
    rows = [0]*n
    for i, (u, v) in enumerate(edge_pairs(n)):
        if mask & (1 << i):
            rows[u] |= 1 << v
            rows[v] |= 1 << u
    return rows


def induced(rows, vertices):
    return sum(((rows[vertices[u]] >> vertices[v]) & 1) << i
               for i, (u, v) in enumerate(edge_pairs(len(vertices))))


def triangle_free(rows):
    return all(not (rows[u] & rows[v]) for u, v in edge_pairs(len(rows))
               if rows[u] & (1 << v))


def representatives():
    labelled = {m for m in range(1 << 15) if triangle_free(adjacency(6, m))}
    count = len(labelled)
    representatives = []
    permutations = list(itertools.permutations(range(6)))
    while labelled:
        first = min(labelled)
        representatives.append(first)
        rows = adjacency(6, first)
        orbit = {induced(rows, order) for order in permutations}
        if not orbit <= labelled:
            raise ValueError("Relabelling orbits are not disjoint")
        labelled.difference_update(orbit)
    return representatives, count


def rootless_key(rows, vertices):
    return min(induced(rows, order) for order in itertools.permutations(vertices))


def rooted_key(rows, vertices):
    a, b, c, d = vertices
    return min(induced(rows, (a, b, c, d)), induced(rows, (a, b, d, c)))


def coordinates():
    raw = {"s0_root0": set(), "s2_root0": set(), "s2_root1": set()}
    roots = set()
    for m in range(8):
        rows = adjacency(3, m)
        if triangle_free(rows):
            raw["s0_root0"].add(rootless_key(rows, range(3)))
    for m in range(64):
        rows = adjacency(4, m)
        if triangle_free(rows):
            raw[f"s2_root{m & 1}"].add(rooted_key(rows, range(4)))
            roots.add(rootless_key(rows, range(4)))
    for root in roots:
        rows = adjacency(4, root)
        raw[f"u4_root{root}"] = {
            profile for profile in range(16)
            if all(not (rows[u] & profile) for u in range(4) if profile & (1 << u))
        }
    return {key: {value: i for i, value in enumerate(sorted(values))}
            for key, values in raw.items()}


def host_counts(mask, indices):
    rows = adjacency(6, mask)
    moments = {key: Counter() for key in indices}
    targets = Counter()
    for order in itertools.permutations(range(6)):
        a, b, c, d, e, f = order
        ab, bc, cd, da, ef = [int(bool(rows[u] & (1 << v)))
                              for u, v in ((a,b), (b,c), (c,d), (d,a), (e,f))]
        targets["cycle"] += ab*bc*cd*da
        targets["matching"] += ab*cd*ef
        targets["edge"] += ab
        block = "s0_root0"
        keys = (rootless_key(rows, (a,b,c)), rootless_key(rows, (d,e,f)))
        moments[block][tuple(indices[block][x] for x in keys)] += 1
        block = f"s2_root{ab}"
        keys = (rooted_key(rows, (a,b,c,d)), rooted_key(rows, (a,b,e,f)))
        moments[block][tuple(indices[block][x] for x in keys)] += 1
        root, reindex = min((induced(rows, tuple(order[i] for i in perm)), perm)
                            for perm in itertools.permutations(range(4)))
        block = f"u4_root{root}"
        keys = [sum(int(bool(rows[leaf] & (1 << order[old]))) << new
                    for new, old in enumerate(reindex)) for leaf in (e,f)]
        moments[block][tuple(indices[block][x] for x in keys)] += 1
    if sum(moments["s0_root0"].values()) != 720:
        raise ValueError("Host probability normalization failed")
    return moments, targets


def verify(certificate):
    epsilon = F(certificate["epsilon"])
    if epsilon != F(1, 50_000_000_000):
        raise ValueError("Unexpected claimed four-cycle inequality")
    indices = coordinates()
    terms = []
    for item in certificate["terms"]:
        block, weight = item["block"], F(item["weight"])
        vector = list(map(F, item["vector"]))
        if weight <= 0 or len(vector) != len(indices[block]):
            raise ValueError("Invalid square weight or dimension")
        terms.append((block, weight, vector))
    hosts, labelled_count = representatives()
    if labelled_count != 5789 or len(hosts) != 38:
        raise ValueError("Unexpected exhaustive graph enumeration")
    slacks = []
    for mask in hosts:
        moments, target = host_counts(mask, indices)
        total = sum(weight*sum(count*vector[i]*vector[j]
                               for (i,j), count in moments[block].items())
                    for block, weight, vector in terms)
        slack = (target["cycle"]-target["matching"]
                 +(F(3,64)+epsilon)*target["edge"]-total)/720
        if slack < 0:
            raise ValueError(f"Negative flag coefficient on host {mask}: {slack}")
        slacks.append({"host_mask": mask, "slack": str(slack)})
    return {"verified": True, "labelled_graphs_examined": 32768,
            "triangle_free_labelled_graphs": labelled_count, "isomorphism_classes": len(hosts),
            "positive_squares": len(terms), "dimensions": {k:len(v) for k,v in indices.items()},
            "host_slacks": slacks,
            "external_certificate_author": "Amir Sarid",
            "external_commit": "ba85f319ebb33ec2194ef4c2929a12916329e1e7"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", default="certificates/sarid-clebsch-tangent.json")
    parser.add_argument("--output", default="results/four-cycle-independent.json")
    args = parser.parse_args()
    raw = Path(args.certificate).read_bytes()
    started = time.monotonic()
    result = verify(json.loads(raw))
    result["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    result["elapsed_seconds"] = time.monotonic()-started
    Path(args.output).write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "host_slacks"}, indent=2))


if __name__ == "__main__":
    main()
