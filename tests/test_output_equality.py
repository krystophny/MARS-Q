"""Exact whole-program output comparison for a performance-only change."""
import os
from pathlib import Path
import unittest
from support import NativeTest, EXE


class OutputEquality(NativeTest):
    @unittest.skipUnless(os.environ.get("CHEASE_REFERENCE_EXECUTABLE"),
                         "set CHEASE_REFERENCE_EXECUTABLE to the parent Make build")
    def test_native_outputs_identical(self):
        reference = Path(os.environ["CHEASE_REFERENCE_EXECUTABLE"])
        before = self.smoke(reference, "reference")
        after = self.smoke(EXE, "candidate")
        self.assertEqual(before.keys(), after.keys())
        print("Comparing exact native bytes:", ", ".join(sorted(before)), flush=True)
        for name in before:
            with self.subTest(file=name):
                self.assertEqual(before[name], after[name], name)
