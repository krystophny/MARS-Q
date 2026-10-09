"""Tiled band factors are bit-identical to independent dense elimination."""
from pathlib import Path
from support import NativeTest, NATIVE, PRECISION, routine


class TiledFactors(NativeTest):
    def test_exact_factors_and_dense_solutions(self):
        files = {"precision.f90": PRECISION, "DECLAR.inc": (NATIVE / "DECLAR.inc").read_text()}
        files["native.f"] = "\n".join(routine(name) for name in ("ALDLT", "LYV", "DWY", "LTXW", "SAXPY", "SDOT"))
        files["oracle.f90"] = Path(__file__).with_name("tiled_factors.f90").read_text()
        self.execute(self.compile(files))
