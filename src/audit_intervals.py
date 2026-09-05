"""Second exact implementation of the five interval certificates.

This module imports no other project module. Generic rational interval
operations evaluate the bounds in Appendix A. The tree is read recursively,
independently of the production checker's stack traversal. Neither checker
formalizes the graph-theoretic argument or establishes research priority.
"""

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("Reversed interval")

    @staticmethod
    def point(x):
        return x if isinstance(x, Interval) else Interval(Q(x), Q(x))

    def __add__(self, other):
        other = self.point(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.point(other))

    def __rsub__(self, other):
        return self.point(other) + (-self)

    def __mul__(self, other):
        other = self.point(other)
        products = [x*y for x in (self.lo, self.hi)
                    for y in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.point(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("Interval division through zero")
        return self * Interval(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return self.point(other) / self

    def square(self):
        if self.lo <= 0 <= self.hi:
            return Interval(Q(0), max(self.lo**2, self.hi**2))
        return Interval(min(self.lo**2, self.hi**2),
                        max(self.lo**2, self.hi**2))

    def minimum(self, other):
        return Interval(min(self.lo, other.lo), min(self.hi, other.hi))


def extension_value(z):
    cap = Q(9, 25)
    if not Q(1, 3) <= z <= Q(1, 2):
        raise ValueError("Invalid extension domain")
    if z <= Q(3, 8):
        return max(Q(9, 128)-z/8,
                   (12*cap-32*z*z-12*z+9)/192)
    return max(z*(Q(1, 2)-z)/2, (cap-8*z*z+3*z)/16)


def leaf_test(box, target, upper_density, cut_cap, lower_density):
    a, mass, degree = box
    c = 1-a
    if mass.lo > (4*cut_cap*c.square()).hi:
        return "empty-mass", Q(0)
    if a.hi <= Q(1, 4):
        extension = (cut_cap-mass/(4*c.square())).hi
        extension += extension_value(1-Q(1, 2)/c.lo)
        if extension <= target:
            return "sparse-set", target-extension
    if (19*degree-18*mass).hi <= 0:
        return "nonpositive-tangent", Q(0)
    baseline = (mass/4).minimum(
        cut_cap+a*(upper_density/2-cut_cap)-(1+a)*mass/(4*c))
    if baseline.hi <= target:
        return "baseline", target-baseline.hi
    if lower_density is not None and a.hi <= Q(1, 4) and mass.lo > 0:
        top = Q(18, 25)*c
        ratio = 2*c-1
        cap = top.minimum(ratio)
        cross = 4*c.square()*(lower_density/2-cut_cap)
        slope = mass*(1-top)-2*mass.square()/top
        if (lower_density >= 2*cut_cap and
                mass.hi <= top.lo**2/2 and cross.hi <= cap.lo and
                slope.lo >= 0):
            curvature = (4*top*cap*mass.square()/
                         (ratio*(cap-cross+2*mass))+
                         slope*cross/ratio-top*mass+2*mass.square())
            # L is increasing on [0, 1/4]. Multiplying the two
            # dependent factors a and 1-a independently would be weaker.
            radius = min(4*a.lo*a.lo, a.lo*(1-a.lo))
            saving = radius*max(curvature.lo, Q(0))/(4*c.hi*c.hi)
            if baseline.hi-saving <= target:
                return "opposite-neighbourhood", target-baseline.hi+saving
    # Only points with 19*d > 18*mu remain. Its positive denominator
    # is at most the interval upper endpoint, even if the lower endpoint
    # is nonpositive. All numerator factors are nonnegative.
    plus = min(c.lo/(1-degree.lo), a.lo/degree.hi)
    minus = min(a.lo/(1-degree.lo), c.lo/degree.hi)
    numerator = (mass*degree.square()*(degree+18*mass)).lo
    denominator_upper = (4*c.square()*(19*degree-18*mass)).hi
    saving = plus*minus*numerator/denominator_upper
    if baseline.hi-saving <= target:
        return "perturbation", target-baseline.hi+saving
    return None, baseline.hi-saving


def replay(tree, upper_density, cut_cap, lower_density=None,
           target=Q(2543, 100000)):
    if upper_density/2 < cut_cap:
        raise ValueError("Density row outside the audited scope")
    cursor = 0
    reasons = Counter()
    margins = {}
    maximum_depth = 0

    def visit(box, depth):
        nonlocal cursor, maximum_depth
        if cursor == len(tree):
            raise ValueError("Missing tree node")
        symbol = tree[cursor]
        cursor += 1
        maximum_depth = max(maximum_depth, depth)
        if symbol == "L":
            reason, margin = leaf_test(box, target, upper_density,
                                       cut_cap, lower_density)
            if reason is None:
                raise ValueError(f"Unproved leaf: {box}, bound {margin}")
            reasons[reason] += 1
            margins[reason] = min(margins.get(reason, margin), margin)
            return
        if symbol not in "012":
            raise ValueError("Invalid tree symbol")
        axis = int(symbol)
        old = box[axis]
        midpoint = (old.lo+old.hi)/2
        for piece in (Interval(old.lo, midpoint), Interval(midpoint, old.hi)):
            child = list(box)
            child[axis] = piece
            visit(tuple(child), depth+1)

    visit((Interval(Q(0), Q(1, 2)), Interval(Q(0), 4*cut_cap),
           Interval(Q(0), Q(1))), 0)
    if cursor != len(tree):
        raise ValueError("Trailing tree nodes")
    return {"nodes": cursor, "leaves": sum(reasons.values()),
            "maximum_depth": maximum_depth, "leaf_reasons": dict(reasons),
            "minimum_margins": {k: str(v) for k, v in margins.items()}}


def audit(root):
    rows = [
        ("base", Q(0), Q(3, 10), Q(2, 47), False),
        ("band1", Q(3, 10), Q(31, 100), Q(2, 47), True),
        ("band2", Q(31, 100), Q(63, 200), Q(2, 47), True),
        ("band3", Q(63, 200), Q(3197, 10000), Q(2, 47), True),
        ("middle", Q(3197, 10000), Q(7, 20), Q(1, 25), False),
    ]
    last = Q(0)
    reports = []
    for name, lo, hi, cut, enabled in rows:
        if lo != last or hi <= lo:
            raise ValueError("Gap in density coverage")
        if cut == Q(1, 25) and lo < Q(3197, 10000):
            raise ValueError("Inapplicable density-specific max-cut theorem")
        last = hi
        data = (root/f"results/opposite-{name}.tree").read_bytes()
        entry = replay(data.decode("ascii").strip(), hi, cut, lo if enabled else None)
        entry.update(certificate=f"opposite-{name}", sha256=hashlib.sha256(data).hexdigest(),
                     density_interval=[str(lo), str(hi)], cut_cap=str(cut))
        reports.append(entry)
    target = Q(2543, 100000)
    independence = Q(9, 128)-Q(9, 25)/8
    high = last/8-last**2/4+(Q(3, 64)+Q(1, 50_000_000_000))/4
    if last != Q(7, 20) or not (independence < target and high < target):
        raise ValueError("Uncovered complementary regime")
    return {"verified": True, "target": str(target), "rows": reports,
            "total_nodes": sum(x["nodes"] for x in reports),
            "total_leaves": sum(x["leaves"] for x in reports),
            "independence_bound": str(independence), "high_density_bound": str(high),
            "scope": "second exact interval implementation; no other project imports",
            "human_mathematical_review": False, "end_to_end_formal_proof": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    report = audit(args.root)
    (args.root/"results/interval-audit.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
