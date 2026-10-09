"""Factorization status, edge inputs and independent dense solutions."""
from pathlib import Path
from support import NativeTest, NATIVE, PRECISION, routine


class BandEdges(NativeTest):
    def test_edge_contract_and_dense_solution(self):
        files = {"precision.f90": PRECISION, "DECLAR.inc": (NATIVE / "DECLAR.inc").read_text()}
        files["native.f"] = "\n".join(routine(name) for name in ("ALDLT", "LYV", "DWY", "LTXW", "SAXPY", "SDOT"))
        files["oracle.f90"] = Path(__file__).with_name("band_edges.f90").read_text()
        self.execute(self.compile(files))
