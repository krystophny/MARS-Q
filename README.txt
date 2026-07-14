MarsQ: source code for MARS-F/K/Q, run "make" under Linux to compile, update makefile for proper compiler
CheaseMerge: source code for CHEASE, solve fixed boundary Grad-Shafranov equation to provide input to MARS-*
EXAMPLE: various test examples
RZplot: post-processing Matlab scripts
Manual: MARS-F/K/Q user manual




Portable local builds
=====================

The `MarsQ_2FK/makefile` provides explicit compiler targets which leave the
historical tracked `marsq.x` untouched and write executables under `build/`:

    make -C MarsQ_2FK gnu
    source /opt/intel/oneapi/setvars.sh
    make -C MarsQ_2FK ifx
    export PATH=/path/to/nvhpc/comm_libs/hpcx/bin:/path/to/nvhpc/compilers/bin:$PATH
    make -C MarsQ_2FK nvhpc

The GNU target includes the compatibility flags required by this legacy source.
The NVIDIA target intentionally uses `-O1`; nvfortran 26.5 crashes internally
while compiling `kinetic.f` at `-O2`.
The equivalent-current (`KKF=-3`) Biot-Savart kernel uses OpenMP; set
`OMP_NUM_THREADS` to the desired physical-core count at runtime.
