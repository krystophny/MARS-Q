"""Regression coverage for the local MARS portability and NTV changes.

The source-contract tests run without a compiler.  Set ``MARS_EXE`` to add
black-box input-validation tests against a built MARS executable.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
MARS_SOURCE = (ROOT / "MarsQ_2FK" / "marsq.f").read_text()
KINETIC_SOURCE = (ROOT / "MarsQ_2FK" / "kinetic.f").read_text()
KINETIC_MODULE = (ROOT / "MarsQ_2FK" / "kineticm.f").read_text()
TORQUE_SOURCE = (ROOT / "MarsQ_2FK" / "torque.f").read_text()
PAMS_SOURCE = (ROOT / "MarsQ_2FK" / "pams.f").read_text()
NEWRUN = (ROOT / "MarsQ_2FK" / "newrun.inc").read_text()
MAKEFILE = (ROOT / "MarsQ_2FK" / "makefile").read_text()
CHEASE_MAKEFILE = (ROOT / "CheaseMerge" / "makefile").read_text()
BUILD_SPEC = importlib.util.spec_from_file_location(
    "build_with_provenance", ROOT / "tools" / "build_with_provenance.py"
)
assert BUILD_SPEC is not None and BUILD_SPEC.loader is not None
build_provenance = importlib.util.module_from_spec(BUILD_SPEC)
BUILD_SPEC.loader.exec_module(build_provenance)


def executable() -> Path | None:
    value = os.environ.get("MARS_EXE")
    if not value:
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise AssertionError(f"MARS_EXE does not exist: {path}")
    return path


def run_with_input(run_input: str) -> subprocess.CompletedProcess[str]:
    exe = executable()
    if exe is None:
        raise unittest.SkipTest("set MARS_EXE for executable input tests")
    with tempfile.TemporaryDirectory(prefix="mars-regression-") as tmp:
        Path(tmp, "RUN.IN").write_text(textwrap.dedent(run_input))
        return subprocess.run(
            [str(exe)],
            cwd=tmp,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=20,
            check=False,
        )


def minimal_input(*, kinetic: str = "", qlin: str = "", outopt: str = "") -> str:
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


class SourceContractTests(unittest.TestCase):
    """Fast tests covering every local change since upstream 8824bb1."""

    def test_portable_compiler_targets_are_out_of_tree(self) -> None:
        for profile in ("gnu", "gnu-debug", "ifx", "nvhpc"):
            self.assertIn(f"{profile}:", MAKEFILE)
            self.assertIn(f"--profile {profile}", MAKEFILE)
        self.assertIn(
            "-fallow-argument-mismatch",
            build_provenance.PROFILES["gnu"]["flags"],
        )
        self.assertIn("-fcheck=all", build_provenance.PROFILES["gnu-debug"]["flags"])
        self.assertIn("-fbacktrace", build_provenance.PROFILES["gnu-debug"]["flags"])
        self.assertIn("-qopenmp", build_provenance.PROFILES["ifx"]["flags"])
        self.assertEqual(
            build_provenance.PROFILES["nvhpc"]["flags"],
            ["-O1", "-r8", "-mp", "-Mextend"],
        )
        self.assertIn(".NOTPARALLEL: gnu gnu-debug ifx nvhpc", MAKEFILE)
        self.assertIn("rm -f *.o *.mod *.d lsode/*.o lsode/*.mod", MAKEFILE)

    def test_build_manifest_binds_the_record_to_the_binary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="mars-manifest-") as temporary:
            directory = Path(temporary)
            binary = directory / "marsq-test.x"
            binary.write_bytes(b"test MARS executable")
            manifest_path = directory / "marsq-test.x.provenance.json"
            manifest = {
                "schema": build_provenance.SCHEMA,
                "artifact": build_provenance.artifact_record(binary),
            }
            build_provenance.write_manifest(manifest_path, manifest)
            validated = build_provenance.validate_manifest(manifest_path)
            self.assertEqual(
                validated["artifact"]["sha256"],
                build_provenance.sha256_file(binary),
            )
            binary.write_bytes(b"tampered MARS executable")
            with self.assertRaisesRegex(ValueError, "mismatch"):
                build_provenance.validate_manifest(manifest_path)

    def test_linked_libraries_use_otool_on_darwin(self) -> None:
        binary = Path("/tmp/marsq-test.x")
        with (
            mock.patch.object(
                build_provenance.platform, "system", return_value="Darwin"
            ),
            mock.patch.object(
                build_provenance.shutil, "which", return_value="/usr/bin/otool"
            ),
            mock.patch.object(
                build_provenance,
                "run_output",
                return_value=(
                    "/tmp/marsq-test.x:\n"
                    "\t/usr/lib/libSystem.B.dylib (compatibility version 1.0.0)"
                ),
            ) as run_output,
        ):
            command, libraries = build_provenance.linked_libraries(binary)

        self.assertEqual(command, ["/usr/bin/otool", "-L", str(binary)])
        self.assertEqual(len(libraries), 2)
        run_output.assert_called_once_with(command, required=False)

    def test_linked_libraries_use_ldd_on_linux(self) -> None:
        binary = Path("/tmp/marsq-test.x")
        with (
            mock.patch.object(
                build_provenance.platform, "system", return_value="Linux"
            ),
            mock.patch.object(
                build_provenance.shutil, "which", return_value="/usr/bin/ldd"
            ),
            mock.patch.object(
                build_provenance, "run_output", return_value="libgfortran.so.5"
            ) as run_output,
        ):
            command, libraries = build_provenance.linked_libraries(binary)

        self.assertEqual(command, ["/usr/bin/ldd", str(binary)])
        self.assertEqual(libraries, ["libgfortran.so.5"])
        run_output.assert_called_once_with(command, required=False)

    def test_linked_libraries_allow_missing_inspection_tool(self) -> None:
        with mock.patch.object(build_provenance.shutil, "which", return_value=None):
            command, libraries = build_provenance.linked_libraries(Path("marsq.x"))

        self.assertIsNone(command)
        self.assertEqual(libraries, [])

    def test_build_profiles_record_compiler_and_exact_flags(self) -> None:
        for profile in build_provenance.PROFILES.values():
            self.assertTrue(profile["binary"].startswith("marsq-"))
            self.assertTrue(profile["compiler"])
            self.assertGreater(len(profile["flags"]), 3)

    def test_mpi_wrapper_path_preserves_argv_zero_dispatch(self) -> None:
        wrapper = shutil.which("mpif90")
        if wrapper is None:
            self.skipTest("mpif90 is not available")
        self.assertEqual(
            build_provenance.resolve_program("mpif90"), Path(wrapper).absolute()
        )

    def test_chease_targets_preserve_required_legacy_initialization(self) -> None:
        self.assertIn("$(BUILD_DIR)/chease-gnu.x", CHEASE_MAKEFILE)
        self.assertIn("-finit-local-zero", CHEASE_MAKEFILE)
        self.assertIn("-mcmodel=medium", CHEASE_MAKEFILE)
        self.assertIn("$(BUILD_DIR)/chease-ifx.x", CHEASE_MAKEFILE)
        self.assertIn("-init=zero", CHEASE_MAKEFILE)

    def test_namelist_reads_report_portable_diagnostics(self) -> None:
        for group in ("BASIC", "FEEDBACK", "KINETIC", "QLIN", "NUMERIC", "OUTOPT"):
            self.assertIn(f"READ(CHNAME,{group},IOSTAT=IOS,IOMSG=IOMSG)", MARS_SOURCE)
            self.assertIn(f"ERROR AT NAMELIST READ: {group}", MARS_SOURCE)
        self.assertIn("IF (IGO.LT.0) STOP 1", MARS_SOURCE)

    def test_input_bnm_header_is_checked_before_integer_conversion(self) -> None:
        self.assertIn("READ(CHOUTP,*) TEMP1,TEMP2", MARS_SOURCE)
        self.assertIn("MAX_EC  = NINT(TEMP1)", MARS_SOURCE)
        self.assertIn("INPUT_BNM HEADER MUST CONTAIN INTEGERS", MARS_SOURCE)
        self.assertIn("IBNM_EC.LT.2.OR.IBNM_EC.GE.NTOT", MARS_SOURCE)

    def test_equivalent_current_kernel_has_complete_openmp_scoping(self) -> None:
        self.assertIn("!$OMP PARALLEL DO DEFAULT(SHARED)", MARS_SOURCE)
        self.assertIn("PRIVATE(MS,KS,J1,J2,J,CTMP1,CTMP2,TEMP1,TEMP2)", MARS_SOURCE)
        self.assertIn("!$OMP& SCHEDULE(DYNAMIC)", MARS_SOURCE)
        self.assertIn("TEMP2 = TEMP2*SQRT(TEMP2)", MARS_SOURCE)
        self.assertIn("!$OMP END PARALLEL DO", MARS_SOURCE)

    def test_mars_k_ntv_rejects_unimplemented_cross_terms(self) -> None:
        self.assertIn("KEYTORQ=2 IS NOT IMPLEMENTED FOR KNTV=20/21", MARS_SOURCE)
        self.assertIn("SUM COMPLEX COIL RESPONSES BEFORE COMPUTING NTV", MARS_SOURCE)
        self.assertIn("KEYTORQ=2 CROSS-COUPLING IS IMPLEMENTED ONLY", NEWRUN)

    def test_perturbative_mars_k_contract_is_explicit(self) -> None:
        self.assertIn("KNTV=20/21 WITH INCKIN>0 REQUIRES ODWKCOM=.TRUE.", MARS_SOURCE)
        self.assertIn("FLUID RESPONSE IS FROZEN", MARS_SOURCE)
        self.assertIn("WITH IPERTURB=1", NEWRUN)

    def test_dwk_component_workspace_is_reused_only_when_consistent(self) -> None:
        """KNTV=21 must not allocate the live KJP master workspace twice."""
        allocator = KINETIC_SOURCE[
            KINETIC_SOURCE.index(
                "SUBROUTINE ALLOCATEDWKCOMPMAT"
            ) : KINETIC_SOURCE.index("END SUBROUTINE ALLOCATEDWKCOMPMAT")
        ]
        self.assertIn("IF (ALLOCATED(INDXDWKC)) THEN", allocator)
        self.assertIn("ANY(INDXDWKC.NE.EXPECTED)", allocator)
        self.assertIn("ANY(SHAPE(VQ3PERPC).NE.SHAPE(VX1PARAC))", allocator)
        self.assertIn("INCOMPLETE DWK COMPONENT WORKSPACE", allocator)
        self.assertIn("INCONSISTENT DWK COMPONENT MAP", allocator)
        self.assertIn("INCONSISTENT DWK COMPONENT WORKSPACE", allocator)
        self.assertLess(
            allocator.index("IF (ALLOCATED(INDXDWKC)) THEN"),
            allocator.index("ALLOCATE (INDXDWKC(NSPECIES,5))"),
        )
        self.assertLess(
            allocator.index("RETURN", allocator.index("IF (ALLOCATED(INDXDWKC))")),
            allocator.index("ALLOCATE (INDXDWKC(NSPECIES,5))"),
        )

    def test_dwk_surface_workspace_is_serialized_in_production(self) -> None:
        """ODWK component files must not depend on optional SFD diagnostics."""
        coefficient = KINETIC_SOURCE[
            KINETIC_SOURCE.index("SUBROUTINE KJPCOEFF") : KINETIC_SOURCE.index(
                "SUBROUTINE KJPFILL"
            )
        ]
        calculator = KINETIC_SOURCE[
            KINETIC_SOURCE.index("SUBROUTINE CALCDWKCOMP") : KINETIC_SOURCE.index(
                "END SUBROUTINE CALCDWKCOMP"
            )
        ]
        self.assertIn("INCSFD = 0", coefficient)
        self.assertIn(
            "IF (ODWKCOM) CALL WRITE_SURFACE_QUANTITIES(JS,KGRID)",
            coefficient,
        )
        self.assertNotIn("IF (INCSFD.GT.0) CALL WRITE_SURFACE_QUANTITIES", coefficient)
        self.assertIn("CALL READ_SURFACE_QUANTITIES (1,2)", calculator)
        self.assertIn("CALL READ_SURFACE_QUANTITIES (IS+1,1)", calculator)
        self.assertIn("DWK CACHE/FIELD/OPPARA/OPPERP/PARA/PERP MAXIMA", calculator)
        self.assertIn("INVALID ZERO DWK CONTRACTION INPUT", calculator)

    def test_frequency_scratch_is_thread_private(self) -> None:
        """KBTIME's per-surface RUU2 must not race across OpenMP workers."""
        self.assertIn("THREADPRIVATE( RCHIHK,RSS,RUU,RUU2 )", KINETIC_MODULE)
        self.assertIn("LAM/RUU2", KINETIC_SOURCE)

    def test_external_frozen_perturbation_import_is_strict_and_atomic(self) -> None:
        self.assertIn("KPERTREAD", NEWRUN)
        self.assertIn("SUBROUTINE READPERTURB", MARS_SOURCE)
        self.assertIn("BPLASMA_INPUT HEADER MISMATCH", MARS_SOURCE)
        self.assertIn("XPLASMA_INPUT HEADER MISMATCH", MARS_SOURCE)
        self.assertIn("BPLASMA_INPUT HAS TRAILING DATA", MARS_SOURCE)
        self.assertIn("XPLASMA_INPUT HAS TRAILING DATA", MARS_SOURCE)
        self.assertIn("BPLASMA_INPUT HAS MALFORMED TRAILING DATA", MARS_SOURCE)
        self.assertIn("XPLASMA_INPUT HAS MALFORMED TRAILING DATA", MARS_SOURCE)
        self.assertIn("BPLASMA_INPUT HAS NON-FINITE DATA", MARS_SOURCE)
        self.assertIn("XPLASMA_INPUT HAS NON-FINITE DATA", MARS_SOURCE)
        self.assertIn("XPLASMA_INPUT PROFILE COLUMNS DISAGREE", MARS_SOURCE)
        self.assertIn("EXTERNAL XPLASMA PROFILE RELATIVE L2", MARS_SOURCE)
        self.assertNotIn("NINT(FN).NE.NINT(RNTOR)", MARS_SOURCE)
        self.assertIn("ASSIGN ONLY AFTER BOTH FILES HAVE PASSED", MARS_SOURCE)
        self.assertLess(
            MARS_SOURCE.index("CALL READPERTURB"),
            MARS_SOURCE.index("CALL OUTPUT(ISWEEP,"),
        )
        self.assertIn("KEYTORQ.EQ.2.OR.KPERTREAD.EQ.1", MARS_SOURCE)
        self.assertEqual(MARS_SOURCE.count("CALL READPERTURB"), 4)
        self.assertIn(
            "IF (INCFEED.GE.0.AND.KPERTREAD.NE.1) CALL FEEDOUT",
            MARS_SOURCE,
        )
        self.assertLess(
            MARS_SOURCE.index("CALL FEEDOUT"),
            MARS_SOURCE.index("CALL TORQNTV", MARS_SOURCE.index("CALL FEEDOUT")),
        )
        self.assertLess(
            MARS_SOURCE.index("CALL TORQNTV", MARS_SOURCE.index("CALL FEEDOUT")),
            MARS_SOURCE.rindex("CALL READPERTURB"),
        )
        self.assertIn("EXTERNAL B/X FIELD", MARS_SOURCE)
        self.assertIn("1172 FORMAT(14(E24.16E3,1X))", MARS_SOURCE)

    def test_keeptfun_is_limited_to_the_kinetic_kjp_window(self) -> None:
        """The carrier fix must not alter the discarded feedback solve."""
        self.assertIn("KEEPTFUN", NEWRUN)
        self.assertIn("KTREST", (ROOT / "MarsQ_2FK" / "feedbackm.f").read_text())
        self.assertIn("KTREST = 1", (ROOT / "MarsQ_2FK" / "feedback.f").read_text())
        window = MARS_SOURCE[
            MARS_SOURCE.index("WRITE (*,'(\"CALLING KJP\")'") : MARS_SOURCE.index(
                "MULTIPLY EQUATIONS INSIDE PLASMA BY EQFAC"
            )
        ]
        self.assertIn("IF (KTREST.NE.0)", window)
        self.assertEqual(window.count("T(JKT)  = TSAVE(JKT)"), 1)
        self.assertEqual(window.count("T(JKT)  = 1.0"), 1)
        self.assertLess(window.index("T(JKT)  = TSAVE(JKT)"), window.index("CALL KJP"))
        self.assertLess(window.index("CALL KJP"), window.index("T(JKT)  = 1.0"))

    def test_kinetic_equilibrium_filter_is_explicitly_selectable(self) -> None:
        """KSMOOTHB exposes the edge-spectrum diagnostic without changing defaults."""
        globalm = (ROOT / "MarsQ_2FK" / "globalm.f").read_text()
        self.assertIn("KSMOOTHB", NEWRUN)
        self.assertIn("KSMOOTHB  = 1", MARS_SOURCE)
        self.assertIn("NKSMOOTHB.LT.0", MARS_SOURCE)
        self.assertIn("KSMOOTHB.NE.0.AND.KSMOOTHB.NE.1", MARS_SOURCE)
        self.assertIn("KSMOOTH = KSMOOTHB", KINETIC_SOURCE)
        self.assertIn("NKSMOOTHB,NKSMOOTHR,NKSINGULAR,KSMOOTHB", globalm)

    def test_equilibrium_drift_factor_is_rebuilt_after_b_filter(self) -> None:
        """The drift derivative must use the same filtered B spectrum as H."""
        marker = "C     SMOOTH EQUILIBRIUM B FIELD FOR BAD EQUILIBRIUM"
        start = KINETIC_SOURCE.index(marker)
        end = KINETIC_SOURCE.index("      DO J=1,NCHI\n         DO JS=1,NRP1", start)
        window = KINETIC_SOURCE[start:end]
        self.assertIn("BPK(JS,J,1)=RJA(JS,J)**2*BK(JS,J,1)", window)
        self.assertIn("BPK(JS,J,2)=RJAM(JS,J)**2*BK(JS,J,2)", window)
        self.assertLess(window.index("ENDIF"), window.index("BPK(JS,J,1)"))

    def test_frozen_field_ntv_skips_feedback_postprocessing(self) -> None:
        """External B/X torque must not enter unrelated feedback diagnostics."""
        output = MARS_SOURCE[
            MARS_SOURCE.index("CYQLIU 15/04/1999") : MARS_SOURCE.index(
                "CYQLIU 12/05/2003"
            )
        ]
        self.assertIn("KPERTREAD.NE.1", output)
        self.assertLess(output.index("CALL FEEDOUT"), output.index("CALL READPERTURB"))
        self.assertLess(output.index("CALL READPERTURB"), output.index("CALL TORQNTV"))

    def test_perturbative_ntv_builds_passive_operator_before_output(self) -> None:
        """Native and imported perturbations need reciprocal DWK blocks."""
        energy_guard = MARS_SOURCE[
            MARS_SOURCE.index(
                "IF (NCASE.NE.6.AND.NCASE.NE.10.AND.KEFORM.NE.0"
            ) : MARS_SOURCE.index("WRITE(*,*) 'AFTER ENERGYMAT'")
        ]
        self.assertIn("KPERTREAD.NE.1", energy_guard)
        self.assertIn("CALL ENERGYMAT", energy_guard)
        pre_output_prepare = MARS_SOURCE[
            MARS_SOURCE.index("IF (KPERTREAD.EQ.1) CALL READPERTURB") :
            MARS_SOURCE.index("CALL OUTPUT(ISWEEP,")
        ]
        self.assertIn(
            "IF (KNTV.EQ.21.AND.INCKIN.GT.0.AND.ISWEEP.EQ.NSWEEP) THEN",
            pre_output_prepare,
        )
        self.assertNotIn("IPERTURB.EQ.1", pre_output_prepare)
        self.assertEqual(
            pre_output_prepare.count("CALL PREPAREKINETICENERGYMAT("), 1
        )
        self.assertGreaterEqual(pre_output_prepare.count("CALL READPERTURB"), 2)
        self.assertLess(
            pre_output_prepare.index("CALL PREPAREKINETICENERGYMAT("),
            pre_output_prepare.rindex("IF (KPERTREAD.EQ.1) CALL READPERTURB"),
        )
        prepare = KINETIC_SOURCE[
            KINETIC_SOURCE.index(
                "SUBROUTINE PREPAREKINETICENERGYMAT"
            ) : KINETIC_SOURCE.index("SUBROUTINE ENERGYMAT")
        ]
        self.assertIn("KJPKEY = 0", prepare)
        self.assertIn("KPBKEY = 1", prepare)
        self.assertIn("CALL LINEAR(ASUBM,BSUBM,CSUBM,DSUBM,", prepare)
        linear = MARS_SOURCE[
            MARS_SOURCE.index("SUBROUTINE LINEAR(") : MARS_SOURCE.index(
                "SUBROUTINE PLASMALIN("
            )
        ]
        self.assertIn("IF (.NOT.(KPERTREAD.EQ.1.AND.KJPKEY.EQ.0))", linear)
        self.assertIn("CALL FEEDM(ASUBM,BSUBM,CSUBM,DSUBM,", linear)
        linear_kjp = MARS_SOURCE[
            MARS_SOURCE.index("WRITE (*,'(\"CALLING KJP\")'") : MARS_SOURCE.index(
                "MULTIPLY EQUATIONS INSIDE PLASMA BY EQFAC"
            )
        ]
        self.assertIn("KDWKREAD.EQ.1.AND.KJPKEY.EQ.0", linear_kjp)
        self.assertIn("CALL KJP(ASUBM,BSUBM,CSUBM,DSUBM,", linear_kjp)
        energy = KINETIC_SOURCE[
            KINETIC_SOURCE.index("SUBROUTINE ENERGYMAT") : KINETIC_SOURCE.index(
                "SUBROUTINE KDWFMAGP",
                KINETIC_SOURCE.index("SUBROUTINE ENERGYMAT"),
            )
        ]
        self.assertNotIn("CALL PREPAREKINETICENERGYMAT(", energy)
        final_energy = MARS_SOURCE[
            MARS_SOURCE.index("IF (NCASE.NE.6.AND.NCASE.NE.10.AND.KEFORM.NE.0") :
            MARS_SOURCE.index("WRITE(*,*) 'AFTER ENERGYMAT'")
        ]
        self.assertIn(
            "IF (.NOT.(KNTV.EQ.21.AND.INCKIN.GT.0))", final_energy
        )
        self.assertLess(
            final_energy.index("CALL PREPAREKINETICENERGYMAT("),
            final_energy.index("CALL ENERGYMAT"),
        )
        self.assertLess(
            energy.index("CALL CALCDWKCOMP("),
            energy.index("IF (KENORM.EQ.2) THEN"),
        )

    def test_mars_k_component_contraction_receives_all_matrices(self) -> None:
        """The legacy implicit interface must never permit a bare DWK call."""
        self.assertNotIn("CALL CALCDWKCOMP\n", KINETIC_SOURCE + TORQUE_SOURCE)

    def test_hot_ion_frequency_diagnostic_guards_species_three(self) -> None:
        """A two-species thermal run has no SLAM0(:,3) diagnostic value."""
        diagnostic = KINETIC_SOURCE[
            KINETIC_SOURCE.index("LAMH = 0.") :
            KINETIC_SOURCE.index("FREQK(JS,1)  = AOMEGABPN")
        ]
        self.assertIn("IF (NSPECIES.GE.3) THEN", diagnostic)
        self.assertLess(
            diagnostic.index("IF (NSPECIES.GE.3) THEN"),
            diagnostic.index("LAMH = SLAM0(L,3)"),
        )
        self.assertIn(
            "CALL CALCDWKCOMP(ASUBM,BSUBM,CSUBM,DSUBM,", TORQUE_SOURCE
        )
        self.assertIn(
            "SUBROUTINE TORQNTV(ASUBM,BSUBM,CSUBM,DSUBM,", TORQUE_SOURCE
        )
        self.assertEqual(MARS_SOURCE.count("CALL TORQNTV(ASUBM"), 2)
        self.assertIn("CALL TORQNTV(A,B,C,D,E,F,G,H)", PAMS_SOURCE)
        self.assertIn("SUBROUTINE OUTPUT(ISW,", MARS_SOURCE)
        self.assertIn("CALL OUTPUT(ISWEEP,", MARS_SOURCE)

    def test_dwk_cache_restart_is_explicit_and_narrow(self) -> None:
        self.assertIn("KDWKREAD", NEWRUN)
        self.assertIn("KDWKREAD MUST BE 0 OR 1", MARS_SOURCE)
        self.assertIn("KDWKREAD=1 REQUIRES KPERTREAD=1, KNTV=21,", MARS_SOURCE)
        self.assertIn("NSWEEP.NE.1", MARS_SOURCE)
        self.assertIn("ISMPIRUN.NE.0", MARS_SOURCE)
        self.assertIn("IF (KDWKREAD.EQ.1) THEN", KINETIC_SOURCE)
        self.assertIn("KJP: REUSING VALIDATED DWK COMPONENT CACHE", KINETIC_SOURCE)
        self.assertIn("KDWKREAD.NE.1) CALL KDWKDENSITY", TORQUE_SOURCE)


class ExecutableInputTests(unittest.TestCase):
    """Black-box tests that stop before equilibrium files are needed."""

    def test_unknown_namelist_key_is_actionable_and_nonzero(self) -> None:
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

    def test_mars_k_ntv_requires_dwk_components(self) -> None:
        result = run_with_input(
            minimal_input(kinetic="INCKIN=1", qlin="KNTV=20, KEYTORQ=0")
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KNTV=20/21 WITH INCKIN>0 REQUIRES ODWKCOM=.TRUE.", result.stdout)

    def test_mars_k_ntv_cross_term_request_fails_loudly(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1",
                qlin="KNTV=21, KEYTORQ=2",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KEYTORQ=2 IS NOT IMPLEMENTED FOR KNTV=20/21", result.stdout)

    def test_valid_perturbative_ntv_announces_frozen_response(self) -> None:
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

    def test_full_two_species_mars_k_models_parse(self) -> None:
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

    def test_external_import_rejects_unknown_mode(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=2",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KPERTREAD MUST BE 0 OR 1", result.stdout)

    def test_external_import_requires_perturbative_mars_k_ntv(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=0, KPERTREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("WITH INCKIN>0 AND IPERTURB>0", result.stdout)

    def test_external_import_is_available_to_shaing_spectrum_diagnostic(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=0, KPERTREAD=1",
                qlin="KNTV=11, KEYTORQ=1",
            )
        )
        self.assertIn("EXTERNAL FROZEN B/X IMPORT ENABLED", result.stdout)
        self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_dwk_cache_restart_rejects_non_frozen_modes(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KDWKREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KDWKREAD=1 REQUIRES KPERTREAD=1", result.stdout)

    def test_dwk_cache_restart_accepts_strict_frozen_mode(self) -> None:
        result = run_with_input(
            minimal_input(
                kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=1, KDWKREAD=1",
                qlin="KNTV=21, KEYTORQ=0",
                outopt="ODWKCOM=.TRUE.",
            )
        )
        self.assertIn("VALIDATED DWK COMPONENT CACHE ENABLED", result.stdout)
        self.assertNotIn("ERROR AT MARS-K NTV INPUT", result.stdout)

    def test_dwk_cache_restart_rejects_multiple_sweeps(self) -> None:
        run_input = minimal_input(
            kinetic="INCKIN=1, IPERTURB=1, KPERTREAD=1, KDWKREAD=1",
            qlin="KNTV=21, KEYTORQ=0",
            outopt="ODWKCOM=.TRUE.",
        ).replace("&BASIC", "&BASIC\n NSWEEP=2")
        result = run_with_input(run_input)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("KDWKREAD=1 REQUIRES KPERTREAD=1", result.stdout)


if __name__ == "__main__":
    unittest.main()
