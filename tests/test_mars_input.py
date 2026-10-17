"""Black-box MARS input validation for the frozen-field NTV keys.

Runs a built MARS executable (MARS_EXE) on minimal RUN.IN files that stop
before any equilibrium input is needed and checks the validation messages
and exit status.  Skipped when MARS_EXE is not set.  Ported from the fork
ExecutableInputTests (PR #5/#24 lineage) onto the build/tests-harness layout.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest


def executable():
    value = os.environ.get("MARS_EXE")
    if not value:
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise AssertionError("MARS_EXE does not exist: %s" % path)
    return path


def run_with_input(run_input):
    exe = executable()
    if exe is None:
        raise unittest.SkipTest("set MARS_EXE for executable input tests")
    with tempfile.TemporaryDirectory(prefix="mars-input-", dir=os.environ.get("TMPDIR")) as tmp:
        Path(tmp, "RUN.IN").write_text(textwrap.dedent(run_input))
        return subprocess.run([str(exe)], cwd=tmp, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=60, check=False)


def minimal_input(*, kinetic="", qlin="", outopt=""):
    return f"""
        &BASIC
        /
        &FEEDBACK
        /
        &KINETIC
        {kinetic}
        /
        &QLIN
        {qlin}
        /
        &NUMERIC
        /
        &OUTOPT
        {outopt}
        /
    """



class ExecutableInputTests(unittest.TestCase):
    """Black-box tests that stop before equilibrium files are needed."""

    def test_dimensional_thermal_authority_accepts_maxwellian_mixture(self):
        for count in (2, 3):
            with self.subTest(species=count):
                deck = minimal_input(kinetic=f"""
                    KPROFTAUTH=2, NSPECIES={count}, ISPECIES_F0=0,0,0,
                    ESPECIES_Z=1,-1,6, ESPECIES_M=2,.00054463,11.9133904496565248
                """).replace("&BASIC", "&BASIC\n NEXPV=1,NPROFIE=4")
                result = run_with_input(deck)
                # Missing geometry is deliberate: the namelist echo occurs
                # only after the complete input validation has succeeded.
                self.assertIn("&BASIC", result.stdout)
                self.assertNotIn("KPROFTAUTH=2 REQUIRES", result.stdout)

    def test_dimensional_thermal_authority_rejects_unsupported_species(self):
        for count, f0, message in (
            (4, "0,0,0,0", "2 OR 3 SPECIES"),
            (3, "0,0,1", "MAXWELLIAN SPECIES"),
        ):
            with self.subTest(species=count, f0=f0):
                deck = minimal_input(kinetic=f"""
                    KPROFTAUTH=2, NSPECIES={count}, ISPECIES_F0={f0}
                """).replace("&BASIC", "&BASIC\n NEXPV=1,NPROFIE=4")
                result = run_with_input(deck)
                self.assertIn("KPROFTAUTH=2 REQUIRES", result.stdout)
                self.assertIn(message, result.stdout)

    def test_unknown_namelist_key_is_actionable_and_nonzero(self):
        result = run_with_input(minimal_input().replace("&BASIC", "&BASIC\n BADKEY=1"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ERROR AT NAMELIST READ: BASIC", result.stdout)
        self.assertIn("IOMSG=", result.stdout)
        diagnostic = result.stdout.lower()
        self.assertTrue(
            "badkey" in diagnostic
            or "entity name is not member of group" in diagnostic
            or ("line " in diagnostic and "position" in diagnostic),
            diagnostic,
        )

    def test_legacy_kntv20_accepts_disabled_dwk_components(self):
        # KNTV=20 contracts the physical pressure in KDWKDENSITY; it does
        # not need the optional per-component DWK cache used by KNTV=21.
        for outopt in ("", "ODWKCOM=.FALSE."):
            with self.subTest(outopt=outopt):
                result = run_with_input(
                    minimal_input(
                        kinetic="INCKIN=1",
                        qlin="KNTV=20, KEYTORQ=0",
                        outopt=outopt,
                    )
                )
                # Deliberately absent equilibrium data stops the run after
                # namelist validation; its echo proves input was accepted.
                self.assertIn("&BASIC", result.stdout)
                self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_kntv21_requires_dwk_components(self):
        result = run_with_input(
            minimal_input(kinetic="INCKIN=1", qlin="KNTV=21, KEYTORQ=0")
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KNTV=21 WITH INCKIN>0 REQUIRES ODWKCOM=.TRUE.", result.stdout)

    def test_mars_k_ntv_cross_term_request_fails_loudly(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1",
                qlin="KNTV=21, KEYTORQ=2",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KEYTORQ=2 IS NOT IMPLEMENTED FOR KNTV=20/21", result.stdout)

    def test_valid_perturbative_ntv_announces_frozen_response(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        # The deliberately minimal fixture later fails on missing equilibrium
        # input, but it must pass NTV validation and identify frozen-field mode.
        self.assertIn("FLUID RESPONSE IS FROZEN", result.stdout)
        self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_full_two_species_mars_k_models_parse(self):
        common = """
            INCKIN=1, ATAU=(1.0e11,0.0), ALTAU=(0.5,0.0),
            ALPHAP=0.5, ALPHAD=1.0, OMEGACI0=50.0,
            IPERTURB=1, KFASTRUN=1, KENORM=2, KEFORM=2,
            NSPECIES=2, IFOWPSI0=1, ISPECIES_F0=0,0,
            ESPECIES_M=2.0,0.00054463, ESPECIES_Z=1.0,-1.0,
            PSPECIES_AP=1.0,1.0, PSPECIES_AT=1.0,1.0,
            PSPECIES_NP=1.0,1.0, PSPECIES_NTB=1.0,1.0,
            PSPECIES_NTD=1.0,1.0, PSPECIES_NDB=1.0,1.0,
            PSPECIES_FOWP=0.0,0.0, PSPECIES_FOWT=0.0,0.0,
            NPROFK=0
        """
        models = (
            "INUTYPE=0, NPROFUI=1, NPROFUE=1",
            "INUTYPE=1, NUMODEL=1, NPROFUI=3, NPROFUE=3",
        )
        for model in models:
            with self.subTest(model=model):
                result = run_with_input(
                    minimal_input(
                        kinetic=f"{common}, {model}, NUEFFIA=1.0, NUEFFEA=1.0",
                        qlin="KNTV=21, KEYTORQ=0",
                        outopt="ODWKCOM=.TRUE.",
                    )
                )
                self.assertNotIn("ERROR AT NAMELIST READ: KINETIC", result.stdout)
                self.assertIn("FLUID RESPONSE IS FROZEN", result.stdout)

    def test_external_import_rejects_unknown_mode(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=2",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KPERTREAD MUST BE 0 OR 1", result.stdout)

    def test_external_import_requires_perturbative_mars_k_ntv(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=0, KPERTREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("WITH INCKIN>0 AND IPERTURB>0", result.stdout)

    def test_external_import_is_available_to_shaing_spectrum_diagnostic(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=0, KPERTREAD=1",
                qlin="KNTV=11, KEYTORQ=1",
            )
        )
        self.assertIn("EXTERNAL FROZEN B/X IMPORT ENABLED", result.stdout)
        self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_dwk_cache_restart_rejects_non_frozen_modes(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KDWKREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KDWKREAD=1 REQUIRES KPERTREAD=1", result.stdout)

    def test_dwk_cache_restart_accepts_strict_frozen_mode(self):
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=1, KDWKREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertIn("VALIDATED DWK COMPONENT CACHE ENABLED", result.stdout)
        self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_dwk_cache_restart_rejects_multiple_sweeps(self):
        run_input = minimal_input(
            kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=1, KDWKREAD=1",
            qlin="KNTV=21, KEYTORQ=0",
            outopt="ODWKCOM=.TRUE.",
        ).replace("&BASIC", "&BASIC\n NSWEEP=2")
        result = run_with_input(run_input)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KDWKREAD=1 REQUIRES KPERTREAD=1", result.stdout)


