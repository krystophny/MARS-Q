MarsQ: source code for MARS-F/K/Q, run "make" under Linux to compile, update makefile for proper compiler
CheaseMerge: source code for CHEASE, solve fixed boundary Grad-Shafranov equation to provide input to MARS-*
EXAMPLE: various test examples
RZplot: post-processing Matlab scripts
Manual: MARS-F/K/Q user manual

The maintained source of this fork is the repository's default branch:
https://github.com/krystophny/MARS-Q
Use a clean checkout of that branch for released builds. Historical stack
branches and their measurements retain provenance, while new production
artifacts must identify the maintained revision and their explicit input
controls. Commit messages retain the original fork and upstream PR links.

For externally maintained physical kinetic backgrounds, KPROFSPEC=1 reads
PROFSPECIES.IN with explicit species identities, densities and temperature
energies. This frozen-background model is distinct from a self-consistent
equilibrium. KPROFTAUTH=1 remains a two-Maxwellian compatibility path and
requires NPROFIE=4,NEXPV=1. The passing-particle magnetic drift is included by
default (KPASSDRIFT=1). Kinetic runs with a passing resonance on paths that
have no drift term (KFASTRUN=1, IFOW>0, finite birth energy) stop unless
KPASSDRIFT=0 states that the original MARS-K model without it is intended.

Portable builds write executables under build/ and leave the historical
tracked executable untouched:

    make -C MarsQ_2FK gnu
    make -C MarsQ_2FK gnu-debug
    make -C MarsQ_2FK ifx
    make -C MarsQ_2FK nvhpc

Configure the selected compiler's MPI wrapper on PATH first. Each MARS target
cleans compiler intermediates, then writes a binary and its .provenance.json
record with compiler flags, Git source state and artifact SHA-256. Serialize
builds that share a source checkout. Verify an artifact without executing it:

    python3 tools/build_with_provenance.py --verify-manifest /path/to/marsq-gnu.x.provenance.json

The default keeps the dependency-free PAMS backend. Enable the certified
global steady-response backend with a GNU profile and the fortsparse prefix:

    make -C MarsQ_2FK gnu PAMS_CSC=1 FORTSPARSE_PREFIX=/path/to/fortsparse

FORTSPARSE_LIBS supplies any additional platform link flags. The provenance
record also binds the fortsparse modules, library and solver helper.
CheaseMerge supports gnu and ifx targets, retaining the legacy source's
zero-initialization and large-memory-model compiler flags.

Run native regression builds and tests on allocated compute nodes. Configure
CTest from the repository root so that both native source directories and
the complete registered test set are covered:

    cmake -S . -B build/tests
    ctest --test-dir build/tests --output-on-failure

Bind CHEASE_EXECUTABLE and MARS_EXE to freshly built native executables.
MARS_PAMS_FORTSPARSE_PREFIX supplies the existing certified sparse prefix for
the real CSC/UMFPACK caller tests; without it that optional dependency suite
reports a skip, which does not certify a sparse production build.


