import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from audit_intervals import Interval, audit, replay
from verify_perturbation_constant import verify


class PaperAuditTests(unittest.TestCase):
    def test_current_trees_agree_between_implementations(self):
        root = Path(__file__).resolve().parents[1]
        report = audit(root)
        self.assertEqual((report["total_nodes"], report["total_leaves"]), (8565, 4285))
        for row in report["rows"]:
            name = row["certificate"]
            lower, upper = map(Q, row["density_interval"])
            cut = Q(row["cut_cap"])
            enabled = lower if "band" in name else None
            tree = (root/f"results/{name}.tree").read_text().strip()
            self.assertEqual(verify(tree, Q(2543, 100000), upper, cut, Q(9, 25), enabled),
                             row["leaf_reasons"])

    def test_current_certificate_corruptions_are_rejected(self):
        root = Path(__file__).resolve().parents[1]
        tree = (root/"results/opposite-band2.tree").read_text().strip()
        for bad in (tree[:-1], tree+"L", "L", "x"+tree[1:]):
            with self.subTest(corruption=bad[:12]):
                with self.assertRaises(ValueError):
                    replay(bad, Q(63, 200), Q(2, 47), Q(31, 100))
        with self.assertRaisesRegex(ValueError, "Unproved leaf"):
            replay(tree, Q(63, 200), Q(2, 47), Q(31, 100), Q(1, 50))

    def test_signed_interval_arithmetic_and_zero_division(self):
        product = Interval(Q(-2), Q(3))*Interval(Q(-5), Q(7))
        self.assertEqual(product, Interval(Q(-15), Q(21)))
        self.assertEqual(Interval(Q(-2), Q(3)).square(), Interval(Q(0), Q(9)))
        with self.assertRaises(ValueError):
            _ = Interval(Q(1), Q(2))/Interval(Q(-1), Q(1))


if __name__ == "__main__":
    unittest.main()
