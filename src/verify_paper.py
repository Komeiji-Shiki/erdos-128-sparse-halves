"""Replay the numerical ingredients of the current 0.02543 draft.

The written graph-theoretic reductions and the quoted max-cut theorem are
separate mathematical obligations; this is not an end-to-end formal proof.
"""

import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path

from verify_four_cycle import verify as verify_cycle
from verify_perturbation_constant import verify as verify_tree
from sparse_set_bound import independent_bound
from audit_intervals import audit as audit_intervals


def main():
    root = Path(__file__).resolve().parents[1]
    target, alpha = F(2543,100000), F(9,25)
    results = []
    started = time.monotonic()
    bands = (("opposite-base", F(0), F(3,10), F(2,47), None),
             ("opposite-band1", F(3,10), F(31,100), F(2,47), F(3,10)),
             ("opposite-band2", F(31,100), F(63,200), F(2,47), F(31,100)),
             ("opposite-band3", F(63,200), F(3197,10000), F(2,47), F(63,200)),
             ("opposite-middle", F(3197,10000), F(7,20), F(1,25), None))
    previous = F(0)
    for name, left, rho, cut, rho_lower in bands:
        if left != previous or rho <= left or (rho_lower is not None and rho_lower > left):
            raise ValueError("Density partition has a gap or invalid lower bound")
        if cut == F(1,25) and left < F(3197,10000):
            raise ValueError("Density-specific cut theorem used outside its range")
        previous = rho
        path = root/"results"/(name+".tree")
        raw = path.read_bytes()
        counts = verify_tree(raw.decode("ascii").strip(), target, rho, cut, alpha, rho_lower)
        results.append({"certificate":name,"rho_lower":str(left),"rho_upper":str(rho),"cut_upper":str(cut),
                        "leaves":counts,"sha256":hashlib.sha256(raw).hexdigest()})
    certificate_path = root/"certificates/sarid-clebsch-tangent.json"
    cycle = verify_cycle(json.loads(certificate_path.read_bytes()))
    high = F(159,6400)+F(1,200_000_000_000)
    independent = independent_bound(alpha)
    if not high < target or not independent < target:
        raise ValueError("Uncovered density or independence regime")
    crosscheck = audit_intervals(root)
    for original, second in zip(results, crosscheck["rows"], strict=True):
        if (original["sha256"] != second["sha256"] or
                original["leaves"] != second["leaf_reasons"]):
            raise ValueError("The two interval implementations disagree")
    result = {"numerical_ingredients_verified":True,"target":str(target),
              "independence_bound":str(independent),"high_density_bound":str(high),
              "interval_checks":results,"four_cycle_hosts":cycle["isomorphism_classes"],
              "four_cycle_positive_squares":cycle["positive_squares"],
              "second_interval_implementation":crosscheck,
              "elapsed_seconds":time.monotonic()-started,
              "not_checked_by_this_program":["written graph-theoretic reductions",
                  "quoted max-cut theorem", "research novelty", "human authorship and review"],
              "original_problem_solved":False}
    (root/"results/paper-verification.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
