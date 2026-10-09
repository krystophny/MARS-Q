"""Small native-kernel and whole-program regression helpers (standard library only)."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import unittest

ROOT = Path(os.environ.get("CHEASE_ROOT") or Path(__file__).resolve().parents[1])
MARS = (ROOT / "CheaseMerge").is_dir()
NATIVE = ROOT / ("CheaseMerge" if MARS else "src-f90")
EXE = Path(os.environ.get("CHEASE_EXECUTABLE") or NATIVE / ("chease.x" if MARS else "chease"))
PRECISION = "module prec_const\ninteger,parameter::rkind=kind(1.d0)\nend module\n"


def routine(name):
    """Read the production routine, extracting a deck from the fixed-form monolith."""
    if not MARS:
        return (NATIVE / (name.lower() + ".f90")).read_text()
    source = (NATIVE / "chease.f").read_text()
    match = re.search(r"^ +(?:SUBROUTINE|(?:[A-Z*0-9]+ +)?FUNCTION) +" + name + r"\b", source, re.M | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    end = source.find("\nC*DECK", match.end())
    return source[match.start():end if end >= 0 else len(source)]


class NativeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", ROOT))
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def compile(self, files, libraries=()):
        """Compile ordered test fixtures and unmodified production sources."""
        paths = []
        for name, content in files.items():
            path = self.work / name
            path.write_text(content)
            if path.suffix.lower() in (".f", ".f90"):
                paths.append(str(path))
        flags = ["-O0", "-fcheck=all", "-ffree-line-length-none", "-ffixed-line-length-none"]
        if MARS:
            flags += ["-std=legacy", "-fdefault-real-8", "-fdefault-double-8", "-fallow-argument-mismatch", "-mcmodel=medium"]
        result = subprocess.run([os.environ.get("FC", "gfortran"), *flags, *paths, *libraries, "-o", "oracle"],
                                cwd=self.work, capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr)
        return self.work / "oracle"

    def execute(self, exe, *args, success=True):
        result = subprocess.run([str(exe), *map(str, args)], cwd=self.work,
                                capture_output=True, text=True, timeout=300)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def smoke(self, exe=EXE, label="native"):
        """Run an existing repository input deck with the Make-built executable."""
        folder = self.work / label
        folder.mkdir()
        if MARS:
            deck = (ROOT / "EXAMPLE/Equi/Rchease_EXAMPLE_base").read_text().split("<<'EOF'\n", 1)[1].split("\nEOF", 1)[0]
            shutil.copyfile(ROOT / "EXAMPLE/Data/EXPEQ_EXAMPLE", folder / "EXPEQ")
        else:
            deck = (ROOT / "WK/TESTCASES/NTCASE2/chease_namelist_tcase2").read_text()
            (folder / "chease_namelist").write_text(deck)
        (folder / "input").write_text(deck + "\n")
        start = time.monotonic()
        with (folder / "input").open() as stream:
            result = subprocess.run([str(exe)], stdin=stream, cwd=folder,
                                    capture_output=True, text=True, timeout=300)
        print(f"{label}: {time.monotonic()-start:.3f} s", flush=True)
        (folder / "run.log").write_text(result.stdout + result.stderr)
        self.assertEqual(result.returncode, 0, (result.stdout + result.stderr)[-2000:])
        self.assertNotRegex(result.stdout, r"(?i)(no convergence|not converged|negative.*tmf|error stop)")
        outputs = {p.name: p.read_bytes() for p in folder.iterdir()
                   if p.name.startswith(("EQDSK", "EXPEQ.OUT", "NJA", "JSOLVER", "NGA", "NVAC", "NDES", "NSAVE", "RZPEEL")) and p.stat().st_size}
        self.assertTrue(outputs, "CHEASE did not produce numerical output")
        for name, data in outputs.items():
            self.assertNotRegex(data, rb"(?i)(?:^|\s)(?:nan|[-+]?inf)(?:\s|$)", name)
        return outputs
