"""Exact whole-program output comparison for a performance-only change."""
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
import unittest

from support import EXE, ROOT, NativeTest

PREFIXES = ("EQDSK", "EXPEQ.OUT", "NJA", "JSOLVER", "NGA", "NVAC", "NDES", "NSAVE", "RZPEEL")


def unlimited_stack():
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


class OutputEquality(NativeTest):
    def run_example(self, exe, label):
        """Run the shipped example deck and return its numerical output files."""
        folder = self.work / label
        folder.mkdir()
        deck = (ROOT / "EXAMPLE/Equi/Rchease_EXAMPLE_base").read_text()
        deck = deck.split("<<'EOF'\n", 1)[1].split("\nEOF", 1)[0]
        shutil.copyfile(ROOT / "EXAMPLE/Data/EXPEQ_EXAMPLE", folder / "EXPEQ")
        (folder / "input").write_text(deck + "\n")
        start = time.monotonic()
        with (folder / "input").open() as stream:
            result = subprocess.run([str(exe)], stdin=stream, cwd=folder, capture_output=True,
                                    text=True, timeout=300, preexec_fn=unlimited_stack,
                                    env={**os.environ, "OMP_NUM_THREADS": "1"})
        print("%s: %.3f s" % (label, time.monotonic() - start), flush=True)
        self.assertEqual(result.returncode, 0, (result.stdout + result.stderr)[-2000:])
        outputs = {p.name: p.read_bytes() for p in folder.iterdir()
                   if p.name.startswith(PREFIXES) and p.stat().st_size}
        self.assertTrue(outputs, "CHEASE did not produce numerical output")
        return outputs

    @unittest.skipUnless(os.environ.get("CHEASE_REFERENCE_EXECUTABLE") and EXE.exists(),
                         "set CHEASE_EXECUTABLE and CHEASE_REFERENCE_EXECUTABLE (parent Make build)")
    def test_native_outputs_identical(self):
        before = self.run_example(Path(os.environ["CHEASE_REFERENCE_EXECUTABLE"]), "reference")
        after = self.run_example(EXE, "candidate")
        self.assertEqual(before.keys(), after.keys())
        print("Comparing exact native bytes:", ", ".join(sorted(before)), flush=True)
        for name in before:
            with self.subTest(file=name):
                self.assertEqual(before[name], after[name], name)
