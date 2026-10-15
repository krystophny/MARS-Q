#!/usr/bin/env python3
"""Check the branches listed in stacks.toml and build the stack branches.

  tools/build_stack.py            check, then build stack/<name> for every stack
  tools/build_stack.py --check    checks only, no branch is written

Checks: every listed branch has merge-base == base and its diff applies to
base alone; branches of a stack with `pairwise = true` also apply pairwise.
A stack may also list `layers = [{ branch = ..., on = ... }]`: a layer is
based on the tip of its `on` branch (an earlier branch or layer of the same
stack) instead of on base; it must have merge-base == tip of `on` and its
diff on..branch must apply to `on` alone.
Build: stack/<name> = base (or the stack it extends) with the commits of each
branch, then of each layer (on..branch), cherry-picked in manifest order, in
a temporary worktree.  Commits already present are dropped.  Any failed
check or conflict stops with a nonzero exit status and writes nothing.
"""
import argparse
import itertools
import os
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def git(*args, cwd=ROOT, env=None, check=True, inp=None, text=True):
    """Run git; text=False keeps stdout as bytes (CRLF sources in diffs)."""
    r = subprocess.run(["git", *args], cwd=cwd, env=env, input=inp,
                       capture_output=True, text=text)
    if not text:
        r.stderr = r.stderr.decode(errors="replace")
    if check and r.returncode:
        sys.exit("git %s failed:\n%s" % (" ".join(args), r.stderr.strip()))
    return r


def rev(ref):
    return git("rev-parse", "--verify", "--quiet", ref + "^{commit}").stdout.strip()


def applies(base, *branches):
    """True when the diffs base..branch apply one after another to base."""
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=str(Path(tmp) / "index"))
        git("read-tree", base, env=env)
        for b in branches:
            # bytes: universal newlines would strip the CRs of CRLF files
            diff = git("diff", "--binary", base, b, text=False).stdout
            if not diff:
                continue
            r = git("apply", "--cached", "--check", "-", env=env, check=False,
                    inp=diff, text=False)
            if r.returncode:
                return False, r.stderr.strip()
            git("apply", "--cached", "-", env=env, inp=diff, text=False)
    return True, ""


def stack_order(cfg, name):
    st = cfg["stacks"][name]
    parent = st.get("extends")
    return (stack_order(cfg, parent) if parent else []) + list(st["branches"])


def layers(cfg, name=None):
    """(branch, on) pairs of one stack, or of all stacks, in manifest order."""
    names = [name] if name else list(cfg["stacks"])
    return [(l["branch"], l["on"]) for n in names
            for l in cfg["stacks"][n].get("layers", [])]


def check(cfg):
    base = rev(cfg["base"])
    listed = set()
    for name in cfg["stacks"]:
        listed.update(stack_order(cfg, name))
    listed.update(cfg.get("superseded", {}))
    listed.update(cfg.get("pending", {}))
    for deps in cfg.get("depends_on", {}).values():
        listed.update(deps)
    layered = dict(layers(cfg))
    listed -= set(layered)
    ok = True
    for b in sorted(listed):
        if not rev(b):
            print("MISSING  %s" % b); ok = False; continue
        mb = git("merge-base", base, b).stdout.strip()
        n = git("rev-list", "--count", "%s..%s" % (base, b)).stdout.strip()
        good, err = applies(base, b)
        state = "ok" if mb == base and good else "FAIL"
        ok &= state == "ok"
        print("%-4s %-45s commits=%s base=%s %s" % (
            state, b, n, "ok" if mb == base else mb[:9], err))
    for name, st in cfg["stacks"].items():
        known = set(stack_order(cfg, name))
        for b, on in layers(cfg, name):
            if on not in known:
                print("FAIL layer %s: `on` %s is not an earlier branch or layer of %s"
                      % (b, on, name)); ok = False; known.add(b); continue
            known.add(b)
            if not rev(b) or not rev(on):
                print("MISSING  %s (layer on %s)" % (b, on)); ok = False; continue
            tip = rev(on)
            mb = git("merge-base", tip, b).stdout.strip()
            n = git("rev-list", "--count", "%s..%s" % (tip, b)).stdout.strip()
            good, err = applies(tip, b)
            state = "ok" if mb == tip and good else "FAIL"
            ok &= state == "ok"
            print("%-4s %-45s commits=%s on=%s %s" % (
                state, b, n, on if mb == tip else "%s (merge-base %s)" % (on, mb[:9]), err))
    for name, st in cfg["stacks"].items():
        if not st.get("pairwise"):
            continue
        for a, b in itertools.combinations(st["branches"], 2):
            good, err = applies(base, a, b)
            if not good:
                ok = False
                print("FAIL pair %s + %s: %s" % (a, b, err))
        print("pairwise %s: %d pairs checked" % (
            name, len(st["branches"]) * (len(st["branches"]) - 1) // 2))
    return ok


def build(cfg, name):
    base = rev(cfg["base"])
    parent = cfg["stacks"][name].get("extends")
    start = rev("stack/" + parent) if parent else base
    with tempfile.TemporaryDirectory() as tmp:
        wt = Path(tmp) / "wt"
        git("worktree", "add", "--detach", str(wt), start)
        try:
            todo = [(b, base) for b in cfg["stacks"][name]["branches"]]
            todo += [(b, rev(on)) for b, on in layers(cfg, name)]
            for b, since in todo:
                commits = git("rev-list", "--reverse", "%s..%s" % (since, b)).stdout.split()
                # A branch may start with commits of another stacked branch
                # (e.g. feat/* on upstream-clean/*).  Skip those by patch id:
                # re-picking an earlier commit conflicts with a later one.
                present = {l[2:].strip() for l in git(
                    "cherry", "HEAD", b, since, cwd=wt).stdout.splitlines()
                    if l.startswith("- ")}
                for c in commits:
                    if c in present:
                        print("  %s %s: already present (patch id), skipped" % (b, c[:9]))
                        continue
                    before = git("rev-parse", "HEAD", cwd=wt).stdout
                    r = git("cherry-pick", "--empty=drop", c, cwd=wt, check=False)
                    if r.returncode:
                        git("cherry-pick", "--abort", cwd=wt, check=False)
                        sys.exit("CONFLICT building stack/%s at %s (%s):\n%s" % (
                            name, b, c[:9], (r.stdout + r.stderr).strip()))
                    if git("rev-parse", "HEAD", cwd=wt).stdout == before:
                        print("  %s %s: already present, dropped" % (b, c[:9]))
            head = git("rev-parse", "HEAD", cwd=wt).stdout.strip()
        finally:
            git("worktree", "remove", "--force", str(wt), check=False)
    git("branch", "-f", "stack/" + name, head)
    print("stack/%s = %s (%s commits on %s)" % (
        name, head[:9], git("rev-list", "--count", "%s..%s" % (base, head)).stdout.strip(),
        cfg["base"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", default=str(ROOT / "stacks.toml"))
    ap.add_argument("--check", action="store_true", help="run the checks only")
    a = ap.parse_args()
    cfg = tomllib.loads(Path(a.manifest).read_text())
    if not check(cfg):
        sys.exit("checks failed; no stack branch written")
    if a.check:
        return
    done = set()

    def visit(n):
        if n in done:
            return
        parent = cfg["stacks"][n].get("extends")
        if parent:
            visit(parent)
        build(cfg, n)
        done.add(n)

    for n in cfg["stacks"]:
        visit(n)


if __name__ == "__main__":
    main()
