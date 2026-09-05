"""Exact interval certificate for the continuous low-density inequality.

Each leaf certifies an entire rational box, not a point in a numerical grid.
The tree splits its parent into two closed boxes with union equal to the parent.
This checks only the scalar inequality. The graph-theoretic reduction is a
separate written proof, and the external max-cut/four-cycle theorems are not
proved by this program.
"""

import argparse
import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path

from sparse_set_bound import independent_bound, phi
from opposite_neighbourhood import interval_saving


def upper_bound(box, target, rho=F(3197,10000), cut=F(2,47), alpha=None, rho_lower=None):
    (al, ah), (ml, mh), (dl, dh) = box
    cl, ch = 1-ah, 1-al
    if ml > 4*cut*ch**2:
        return "empty-mass", target
    if alpha is not None and ah <= F(1,4):
        zlo = 1-1/(2*cl)
        extension = cut-ml/(4*ch**2)+phi(zlo, alpha)
        if extension <= target:
            return "sparse-set", extension
    if 19*dh <= 18*ml:
        return "nonpositive-tangent", target
    baseline = min(mh/4, cut + ah*(rho/2-cut) - (1+al)*ml/(4*ch))
    if baseline <= target:
        return "baseline", baseline
    if rho_lower is not None and alpha is not None:
        opposite = baseline-interval_saving(box[0],box[1],rho_lower,cut,alpha)
        if opposite <= target:
            return "opposite-neighbourhood", opposite
    # The endpoint radii are positive on the nontrivial interior. Endpoint
    # zeros produce a valid, possibly weak, lower bound of zero.
    plus = min(cl/(1-dl), al/dh)
    minus = min(al/(1-dl), cl/dh)
    saving = plus*minus*ml*dl**2*(dl+18*ml)/(19*dh-18*ml)
    upper = baseline - saving/(4*ch**2)
    return ("perturbation" if upper <= target else "split"), upper


def children(box, axis):
    left, right = list(box), list(box)
    lo, hi = box[axis]
    mid = (lo+hi)/2
    left[axis] = (lo, mid)
    right[axis] = (mid, hi)
    return tuple(left), tuple(right)


def generate(target, max_nodes, rho=F(3197,10000), cut=F(2,47), alpha=None, rho_lower=None):
    root = ((F(0), F(1, 2)), (F(0), 4*cut), (F(0), F(1)))
    stack = [(root, 0)]
    tree = []
    counts = {}
    nodes = max_depth = 0
    started = time.monotonic()
    while stack:
        box, depth = stack.pop()
        nodes += 1
        if nodes > max_nodes:
            raise RuntimeError(f"Node budget exhausted: {len(stack)} boxes unresolved")
        reason, upper = upper_bound(box, target, rho, cut, alpha, rho_lower)
        max_depth = max(max_depth, depth)
        if reason != "split":
            counts[reason] = counts.get(reason, 0)+1
            tree.append("L")
        else:
            # A numerical heuristic may pick the split; acceptance uses only
            # the exact rational bounds above. Every split preserves coverage.
            scales = (0.2, 0.65, 0.035)
            if rho_lower is not None:
                point = tuple(((lo+hi)/2, (lo+hi)/2) for lo,hi in box)
                center_reason, _ = upper_bound(point,target,rho,cut,alpha,rho_lower)
                if center_reason in {"opposite-neighbourhood", "sparse-set", "baseline", "empty-mass"}:
                    # Those candidates do not depend on d. Refining d would
                    # duplicate the same two-dimensional calculation.
                    scales = (0.2, 0.65, 0)
            axis = max(range(3), key=lambda i: scales[i]*float(box[i][1]-box[i][0]))
            tree.append(str(axis))
            left, right = children(box, axis)
            stack.append((right, depth+1))
            stack.append((left, depth+1))
        if nodes % 50000 == 0:
            print(json.dumps({"nodes": nodes, "pending": len(stack),
                              "elapsed_seconds": round(time.monotonic()-started, 3)}), flush=True)
    return "".join(tree), counts, max_depth, time.monotonic()-started


def verify(tree, target, rho=F(3197,10000), cut=F(2,47), alpha=None, rho_lower=None):
    if alpha is not None and (not F(1,3) <= alpha <= F(3,8)
                              or independent_bound(alpha) > target):
        raise ValueError("Independence case is not covered by the target")
    root = ((F(0), F(1, 2)), (F(0), 4*cut), (F(0), F(1)))
    stack = [root]
    counts = {}
    for symbol in tree:
        if not stack:
            raise ValueError("Trailing certificate nodes")
        box = stack.pop()
        if symbol == "L":
            reason, upper = upper_bound(box, target, rho, cut, alpha, rho_lower)
            if reason == "split":
                raise ValueError(f"Unproved leaf: {box}, {upper}")
            counts[reason] = counts.get(reason, 0)+1
        elif symbol in "012":
            left, right = children(box, int(symbol))
            stack.extend((right, left))
        else:
            raise ValueError("Invalid tree symbol")
    if stack:
        raise ValueError("Missing certificate nodes; domain not covered")
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="163/6250")  # 0.02608
    parser.add_argument("--output", default="results/perturbation-interval")
    parser.add_argument("--max-nodes", type=int, default=2_000_000)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--rho", default="3197/10000")
    parser.add_argument("--cut", default="2/47")
    parser.add_argument("--alpha", help="Enable the sparse-set bound with this independence cap")
    parser.add_argument("--rho-lower", help="Enable opposite-side anchors above this density")
    args = parser.parse_args()
    target = F(args.target)
    rho, cut = F(args.rho), F(args.cut)
    alpha = F(args.alpha) if args.alpha else None
    rho_lower = F(args.rho_lower) if args.rho_lower else None
    if alpha is not None and (not F(1,3) <= alpha <= F(3,8)
                              or independent_bound(alpha) > target):
        raise ValueError("Independence case is not covered by the target")
    prefix = Path(args.output)
    tree_path = prefix.with_suffix(".tree")
    if args.verify:
        started = time.monotonic()
        tree = tree_path.read_text(encoding="ascii").strip()
        counts = verify(tree, target, rho, cut, alpha, rho_lower)
        print(json.dumps({"verified": True, "target": str(target), "leaves": counts,
                          "elapsed_seconds": time.monotonic()-started}, indent=2))
        return
    tree, counts, depth, seconds = generate(target, args.max_nodes, rho, cut, alpha, rho_lower)
    # Canonical bytes keep hashes identical on Windows and Unix.
    tree_path.write_bytes((tree+"\n").encode("ascii"))
    result = {"target": str(target), "rho_upper": str(rho), "cut_upper": str(cut),
              "independence_cap": str(alpha) if alpha is not None else None,
              "rho_lower_for_opposite_neighbourhood": str(rho_lower) if rho_lower is not None else None,
              "tangent_parameter": "19/20", "nodes": len(tree), "leaves": counts,
              "maximum_depth": depth, "elapsed_seconds": seconds,
              "tree_sha256": hashlib.sha256(tree_path.read_bytes()).hexdigest(),
              "scope": "continuous scalar inequality; graph reduction supplied separately"}
    prefix.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
