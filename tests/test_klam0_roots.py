"""Manufactured behavioral oracle for the opt-in KLAM0 root policy.

The production Fortran routine receives tabulated frequencies and has one
scalar SLAM0 consumer.  These tests deliberately model that contract: every
strict sign-changing bracket is enumerated, exact nodes/endpoints are not
selected, and the first valid interior root is retained.  The legacy helper
models its first-crossing and ``J.GT.2 .AND. J.LT.NN-1`` guard independently.
"""

from __future__ import annotations

import math
import unittest


def legacy_first(values: list[float], target: float = 0.0) -> int | None:
    """Return the accepted 1-based legacy bracket, if any."""
    n = len(values)
    j = 1
    while j < n:
        if (values[j - 1] - target) * (values[j] - target) > 0.0:
            j += 1
            if j < n:
                continue
        break
    if j > 2 and j < n - 1:
        return j
    return None


def safe_roots(
    values: list[float], target: float = 0.0
) -> tuple[list[tuple[int, float]], tuple[int, float] | None]:
    """Enumerate sign brackets and linearly solve the manufactured spline."""
    roots: list[tuple[int, float]] = []
    for j, (left, right) in enumerate(zip(values, values[1:]), start=1):
        scale = max(1.0, abs(target), abs(left - target), abs(right - target))
        tol = 1.0e-12 * scale
        fleft = left - target
        fright = right - target
        if abs(fleft) <= tol or abs(fright) <= tol:
            continue
        if fleft * fright >= 0.0:
            continue
        fraction = -fleft / (fright - fleft)
        root = float(j - 1) + fraction
        if not 0.0 < root < float(len(values) - 1):
            continue
        roots.append((j, root))
    return roots, roots[0] if roots else None


class Klam0RootPolicyTests(unittest.TestCase):
    def test_safe_path_recovers_first_crossing_rejected_by_legacy_guard(self) -> None:
        values = [1.0, -1.0, 1.0, 1.0, 1.0, 1.0]
        self.assertIsNone(legacy_first(values))
        roots, selected = safe_roots(values)
        self.assertEqual([j for j, _ in roots], [1, 2])
        self.assertIsNotNone(selected)
        assert selected is not None
        self.assertEqual(selected[0], 1)
        self.assertTrue(math.isclose(selected[1], 0.5))

    def test_safe_path_enumerates_multiple_brackets_but_keeps_scalar_contract(self) -> None:
        values = [1.0, -1.0, 1.0, -1.0, 1.0]
        roots, selected = safe_roots(values)
        self.assertEqual([j for j, _ in roots], [1, 2, 3, 4])
        self.assertIsNotNone(selected)
        assert selected is not None
        self.assertEqual(selected[0], 1)

    def test_exact_nodes_and_global_endpoints_are_not_roots(self) -> None:
        values = [0.0, -1.0, 1.0, -1.0, 0.0]
        roots, selected = safe_roots(values)
        self.assertEqual([j for j, _ in roots], [2, 3])
        self.assertIsNotNone(selected)
        assert selected is not None
        self.assertEqual(selected[0], 2)
        self.assertTrue(all(0.0 < root < len(values) - 1 for _, root in roots))

    def test_no_crossing_returns_no_scalar_root(self) -> None:
        roots, selected = safe_roots([1.0, 2.0, 3.0, 4.0])
        self.assertEqual(roots, [])
        self.assertIsNone(selected)


if __name__ == "__main__":
    unittest.main()
