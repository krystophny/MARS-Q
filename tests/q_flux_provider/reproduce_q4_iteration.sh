#!/usr/bin/env bash
# Self-contained review assembly and generated native q4 reproduction.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo 'Usage: reproduce_q4_iteration.sh NEW_OUTPUT_DIR [--prepare-only]' >&2
  exit 2
fi
repro_root=$(realpath -m "$1")
repro_mode=${2:-run}
if [[ "$repro_mode" != run && "$repro_mode" != --prepare-only ]]; then exit 2; fi
if [[ -e "$repro_root" ]]; then echo 'Use a new output directory.' >&2; exit 2; fi
repo_root=$(git rev-parse --show-toplevel)
script_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$repro_root"
cp "$script_root/generate_inputs.py" "$repro_root/generate_inputs.py"
git -C "$repo_root" worktree add --detach "$repro_root/source" 8824bb18e1514b4a27f6357d5fe859d4fe690542 > "$repro_root/worktree.log" 2>&1
# Source starts from the unchanged parent, excluding the unrun provider PR.
for pr in 30 39 41 42; do
  curl --fail --location --silent --show-error "https://github.com/krystophny/MARS-Q/pull/${pr}.diff" -o "$repro_root/pr${pr}.diff"
  git -C "$repro_root/source" apply --check --include=CheaseMerge/chease.f "$repro_root/pr${pr}.diff"
  git -C "$repro_root/source" apply --include=CheaseMerge/chease.f "$repro_root/pr${pr}.diff"
done
sha256sum "$repro_root"/pr*.diff "$repro_root/source/CheaseMerge/chease.f" > "$repro_root/source_seals.txt"
python3 "$repro_root/generate_inputs.py" "$repro_root/cases" --variant mars
make -B -C "$repro_root/source/CheaseMerge" -j2 F95=gfortran \
  F95FLAGS='-O3 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch -mcmodel=medium' \
  LDFLAGS=-mcmodel=medium > "$repro_root/build.log" 2>&1
sha256sum "$repro_root/source/CheaseMerge/chease.x" "$repro_root/cases"/*/EXPEQ "$repro_root/cases"/*/chease_namelist > "$repro_root/run_seals.txt"
if [[ "$repro_mode" == --prepare-only ]]; then
  echo 'Source patches, native build and generated inputs prepared; no native run.'
  exit 0
fi
repro_status=0
for case_name in circular shaped; do
  case_root="$repro_root/cases/$case_name"
  native_status=0
  (cd "$case_root" && env -u CHEASE_Q_OBSERVER_FILE -u CHEASE_GS_OBSERVER_DIR \
    -u CHEASE_NONLINEAR_OBSERVER_DIR OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 timeout 180 "$repro_root/source/CheaseMerge/chease.x" \
    < chease_namelist > chease.log 2>&1) || native_status=$?
  printf '%s\n' "$native_status" > "$case_root/process_returncode.txt"
  python3 - "$case_root" <<'PY' || repro_status=2
import json, math, pathlib, struct, sys
root=pathlib.Path(sys.argv[1]);p=root/'NOUT';rows=[];reason=[]
if not p.exists() or p.stat().st_size==0:
    reason.append('Native NOUT is absent or empty')
else:
    with p.open('rb')as f:
        while h:=f.read(4):
            if len(h)!=4:raise ValueError('truncated record marker')
            n=struct.unpack('<i',h)[0]
            if n<0:raise ValueError('unexpected continued record in bounded case')
            data=f.read(n);tail=f.read(4)
            if len(data)!=n or tail!=h:raise ValueError('invalid native record')
            rows.append(data)
    ns,nt=struct.unpack('<2i',rows[1][:8])
    if(ns,nt)!=(32,32):reason.append('Requested final NS32/NT32 is absent')
    psi=[v[0]for v in struct.iter_unpack('<4d',rows[8])]
    if not psi or max(psi)<=min(psi):
        reason.append('Native poloidal flux has no finite nonzero span')
    for data in rows[2:]:
        if len(data)%8==0:
            if not all(math.isfinite(v[0])for v in struct.iter_unpack('<d',data)):
                reason.append('Native floating payload is nonfinite');break
result={'native_process_returncode':int((root/'process_returncode.txt').read_text()),
        'native_final_grid_nonempty_finite':not reason,'failure_reasons':reason}
(root/'native_disposition.json').write_text(json.dumps(result,indent=2)+'\n')
print(root.name,json.dumps(result))
sys.exit(2 if reason else 0)
PY
done
exit "$repro_status"
