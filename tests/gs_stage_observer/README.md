# Original Grad–Shafranov stages

Opt-in CHEASE_GS_OBSERVER_DIR captures original assembled A before ALDLT,
original B before DIRECT, reduced CPSILI before CENTER, physical CPSICL before
and after relaxation, after MAGAXE, and immediately pre/post SMOOTH. CHECK's
reassembled A/B are separate files. Native arithmetic and data are unchanged;
all files use STATUS=NEW in a precreated unique external directory.

Each matrix generation has its chart and active band. Each nonlinear step has
its original RHS and tagged states; filenames retain both counters. The first
matrix record is int32(NS,NT,N4NSNT,NBAND), then float64 native active band in
Fortran order, excluding NPBAND padding. RHS records have three dimensions and
values; vector records additionally retain the native NUPDWN map. Chart records
contain dimensions(NS,NT,NBPS), sigma/theta and boundary samples/second derivatives.
Metadata contains generation/step, then R0,RZ0,RMAG,RZMAG,SPSIM,R0EXP,B0EXP,RELAX.
A detected-axis move is not assumed to be a chart rebuild.

Registered CTest verifies native stream headers, active/padded original values,
RHS/permutation/state bytes and no caller-matrix mutation, plus the disabled
route. Full enabled-versus-disabled equilibrium bytes must pass separately.
These observations do not establish a physical continuum bound or repair.
