"""`bio`: the BioLang toolchain on its own.

    bio check   FILE [--context k=v]     compile, resolve references, report evidence and confidence
    bio compile FILE [-o out.json]       BioLang → BioIR JSON
    bio run     FILE [--hours H] ...     run the module on the network runtime (or SBML / .bnet)
    bio test    PATH...                  run every .bio file, evaluate its `# test:` lines
    bio repl                             type BioLang, run it, inspect it

`test` lines live in comments so a module stays a plain program:

    # test: TetR final > 5          the last value of a species
    # test: TetR peaks >= 2         number of peaks over the run
    # test: rules >= 3              module facts: rules, entities, unknowns
    # test: confidence >= 0.6       mean confidence over all rules
    # test: alive == 961             organism programs: cells (born), alive, deaths at the last stage;
                                    the program's own `assert:` lines are checked too

GenomeOS is the first application of this toolchain; nothing here needs
the genome layer, so it can be packaged apart (docs/ARCHITECTURE.md).
"""

from __future__ import annotations

import argparse
import math
import operator
import re
import sys
from pathlib import Path
from typing import Any

from genomeos import __version__
from genomeos.ir import Module
from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError

OPS = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
}
TEST_LINE = re.compile(r"^\s*#\s*test:\s*(\S+)\s+(final|peaks|min|max)?\s*(>=|<=|==|!=|>|<)\s*(-?[\d.]+)\s*$")


# ---- test --------------------------------------------------------------------------


def parse_tests(text: str) -> list[tuple[str, str, str, float]]:
    """(subject, measure, op, value) for every `# test:` line."""
    out = []
    for line in text.splitlines():
        m = TEST_LINE.match(line)
        if m:
            out.append((m.group(1), m.group(2) or "final", m.group(3), float(m.group(4))))
    return out


def evaluate(
    module: Module, tests: list[tuple[str, str, str, float]], hours: float = 48.0
) -> list[dict[str, Any]]:
    """Run the module once (if any test needs a trajectory) and judge each test line."""
    from genomeos.runtime import NetworkRuntime

    facts = ("rules", "entities", "unknowns", "confidence", "events")
    body_subjects = ("cells", "alive", "deaths")
    needs_run = any(s not in facts and s not in body_subjects for s, *_ in tests)
    traj = None
    if needs_run and module.rules:
        vm = NetworkRuntime(module, context={}, seed=0)
        traj = vm.run(hours=hours, dt=0.05, initial={})
    results = []
    # an organism program (stages declared): grow it to the last stage; its own `assert:` lines are claims
    body = None
    stages = getattr(module, "stages", None) or []
    if stages and getattr(module, "organism", None) is not None:
        try:
            from genomeos.ir.model import to_minutes
            from genomeos.runtime.body import Body

            horizon = max(to_minutes(st.start, st.unit) for st in stages)
            body = Body(module, seed=0).run(until=horizon)
            for a in body.check_asserts():
                results.append(
                    {"test": f"assert: {a.get('assert')}", "got": a.get("value"), "ok": bool(a.get("ok"))}
                )
        except Exception as e:  # noqa: BLE001  (the organism layer reports its own failure as a claim)
            results.append({"test": "organism grows to its last stage", "got": str(e)[:80], "ok": False})
    for subject, measure, op, value in tests:
        got: float | None
        if subject == "rules":
            got = len(module.rules)
        elif subject == "entities":
            got = len(module.entities)
        elif subject == "unknowns":
            got = len(module.unknowns())
        elif subject == "events":
            got = len(module.events)
        elif subject == "confidence":
            got = sum(r.confidence for r in module.rules) / len(module.rules) if module.rules else 0.0
        elif subject in body_subjects and body is not None:
            summ = body.summary()
            got = {"cells": summ["cells_born"], "alive": summ["alive"], "deaths": summ["deaths"]}[subject]
        elif traj is not None and subject in traj.levels:
            xs = traj.levels[subject]
            got = {"final": xs[-1], "peaks": traj.peaks(subject), "min": min(xs), "max": max(xs)}[measure]
        else:
            got = None
        ok = got is not None and not (isinstance(got, float) and math.isnan(got)) and OPS[op](got, value)
        fact = subject in facts or subject in body_subjects
        label = f"{subject} {op} {value:g}" if fact else f"{subject} {measure} {op} {value:g}"
        results.append({"test": label, "got": got, "ok": ok})
    if traj is not None:
        bad = [s for s in traj.species if any(math.isnan(x) or math.isinf(x) for x in traj.levels[s])]
        results.append({"test": "no NaN or infinite level", "got": bad or "clean", "ok": not bad})
    return results


def cmd_test(args: argparse.Namespace) -> int:
    files: list[Path] = []
    for p in args.paths:
        pp = Path(p)
        files += sorted(pp.rglob("*.bio")) if pp.is_dir() else [pp]
    failed = 0
    total = 0
    for f in files:
        text = f.read_text()
        try:
            module = parse_file(f)
        except BioLangError as e:
            print(f"✗ {f}: does not compile: {e}")
            failed += 1
            total += 1
            continue
        tests = parse_tests(text)
        organism = bool(getattr(module, "stages", None)) and getattr(module, "organism", None) is not None
        results = evaluate(module, tests, hours=args.hours) if (tests or module.rules or organism) else []
        n_ok = sum(1 for r in results if r["ok"])
        total += max(1, len(results))
        failed += sum(1 for r in results if not r["ok"])
        mark = "✓" if all(r["ok"] for r in results) else "✗"
        print(
            f"{mark} {f}: module {module.name}, {len(module.entities)} entities, {len(module.rules)} rules, "
            f"{n_ok}/{len(results)} checks"
        )
        for r in results:
            if not r["ok"] or args.verbose:
                g = r["got"]
                shown = f"{g:.3g}" if isinstance(g, float) else g
                print(f"    {'ok ' if r['ok'] else 'FAIL'} {r['test']}  (got {shown})")
    print(f"{len(files)} files, {total - failed}/{total} checks passed")
    return 1 if failed else 0


# ---- repl --------------------------------------------------------------------------


HELP = """BioLang REPL. Type blocks (gene, protein, rule, ...); a block ends at its closing brace.
  :run [HOURS] [K=V ...]   run the module (default 24 h) and show final levels
  :check                   entities, rules, confidence report
  :show                    the BioIR as JSON
  :load FILE               read a .bio file into the session
  :clear                   start over
  :quit                    leave"""


class Repl:
    def __init__(self, out=sys.stdout) -> None:
        self.lines: list[str] = ["module repl"]
        self.pending: list[str] = []  # blocks that refer to something not declared yet
        self.out = out

    def module(self) -> Module:
        return parse("\n".join(self.lines))

    def say(self, text: str) -> None:
        print(text, file=self.out)

    def handle(self, line: str) -> bool:
        """One input line; returns False to quit."""
        s = line.strip()
        if not s:
            return True
        if s in (":quit", ":q", "exit"):
            return False
        if s in (":help", "?"):
            self.say(HELP)
            return True
        if s == ":clear":
            self.lines = ["module repl"]
            self.say("cleared")
            return True
        if s.startswith(":load "):
            path = Path(s[6:].strip())
            body = [x for x in path.read_text().splitlines() if not x.strip().startswith("module ")]
            self.lines += body
            self.say(f"loaded {path}: {len(body)} lines")
            return self._report()
        if s == ":show":
            import json

            self.say(json.dumps(self.module().to_dict(), indent=1)[:4000])
            return True
        if s == ":check":
            return self._report(verbose=True)
        if s.startswith(":run"):
            return self._run(s[4:].split())
        # a BioLang statement or block; accept it if the module still compiles. A block that
        # refers to something not declared yet waits until the declaration arrives.
        candidate = [*self.lines, s]
        try:
            parse("\n".join(candidate))
        except BioLangError as e:
            if "undeclared" in str(e):
                self.pending.append(s)
                self.say(f"  waiting: {e}")
                return True
            self.say(f"  rejected: {e}")
            return True
        self.lines = candidate
        self._fold_pending()
        return self._report()

    def _fold_pending(self) -> None:
        changed = True
        while changed and self.pending:
            changed = False
            for line in list(self.pending):
                try:
                    parse("\n".join([*self.lines, line]))
                except BioLangError:
                    continue
                self.lines.append(line)
                self.pending.remove(line)
                changed = True

    def _report(self, verbose: bool = False) -> bool:
        m = self.module()
        waiting = f", {len(self.pending)} waiting for a declaration" if self.pending else ""
        self.say(f"  {len(m.entities)} entities, {len(m.rules)} rules, {len(m.events)} events{waiting}")
        if verbose:
            for kind, val in sorted(m.confidence_report().items()):
                self.say(f"    {kind:<10} confidence {val:.2f}")
            for r in m.rules:
                self.say(f"    {r.id}  [{r.evidence.kind.value}] {r.confidence}")
        return True

    def _run(self, argv: list[str]) -> bool:
        from genomeos.runtime import NetworkRuntime

        hours = 24.0
        initial: dict[str, float] = {}
        for a in argv:
            if "=" in a:
                k, v = a.split("=", 1)
                initial[k] = float(v)
            else:
                hours = float(a)
        m = self.module()
        if not m.rules:
            self.say("  nothing to run: no rules yet")
            return True
        vm = NetworkRuntime(m, context={}, seed=0)
        traj = vm.run(hours=hours, dt=0.05, initial=initial)
        for sp in traj.species:
            xs = traj.levels[sp]
            self.say(
                f"  {sp:<14} final {xs[-1]:8.2f}   min {min(xs):8.2f}   max {max(xs):8.2f}   "
                f"peaks {traj.peaks(sp)}"
            )
        return True


def cmd_repl(args: argparse.Namespace) -> int:
    repl = Repl()
    interactive = sys.stdin.isatty()
    if interactive:
        print(f"bio {__version__} · BioLang REPL · :help for commands")
    while True:
        try:
            line = input("bio> " if interactive else "")
        except EOFError:
            break
        if not repl.handle(line):
            break
    return 0


# ---- entry point ---------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    from genomeos.cli import cmd_check, cmd_compile, cmd_run

    ap = argparse.ArgumentParser(
        prog="bio", description="the BioLang toolchain: check, compile, run, test, repl"
    )
    ap.add_argument("--version", action="version", version=f"bio {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check", help="compile and report evidence and confidence")
    p.add_argument("module")
    p.add_argument("--context", nargs="*")
    p.set_defaults(fn=cmd_check)
    p = sub.add_parser("compile", help="BioLang → BioIR JSON")
    p.add_argument("module")
    p.add_argument("-o", "--output")
    p.set_defaults(fn=cmd_compile)
    p = sub.add_parser(
        "run", help="run a module (.bio, .json), an SBML model (.xml) or a Boolean network (.bnet)"
    )
    p.add_argument("module")
    p.add_argument("--hours", type=float, default=48.0)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--context", nargs="*")
    p.add_argument("--init", nargs="*")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--csv")
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("test", help="run every .bio file and its `# test:` lines")
    p.add_argument("paths", nargs="+")
    p.add_argument("--hours", type=float, default=48.0)
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(fn=cmd_test)
    p = sub.add_parser("repl", help="interactive BioLang")
    p.set_defaults(fn=cmd_repl)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
