# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CLI parses the same way on every CPython the project supports, not just on this machine's.

`tests/test_roadmap_work.py::test_the_board_root_can_be_named_by_the_environment` PASSED here and
FAILED on CI at 2355bad (2026-10-02, ubuntu-24.04, `python-version: 3.12`) with

    genomeos: error: unrecognized arguments: a task

It was first read as an argument-order defect in `genomeos work`, and it is not: both orders parse on
this machine. It was then read as a full-suite-only failure caused by lanes writing the shared board
while the test read it, and it is not that either -- the failing line is `assert r.returncode == 0`
on the very first invocation, before any board is read, and the invocation names its own root.

The cause is argparse, and which argparse. `work` declares two optional positionals, `action` and
`task`, and an option may sit between them. `_match_arguments_partial` used to consume BOTH positionals
as soon as it matched one argument string, giving `task` zero of them, so the task string was left over
and the top-level parser rejected it; CPython later added a guard -- where the next argument string is
an option, a trailing positional stays pending rather than being consumed with nothing.
`requires-python` is >=3.12, which does not fix the micro version, so the project does not get to
decide which of the two behaviours its users and CI run.

So the behaviour is carried by `cli.StableParser` instead of being inherited from the interpreter, and
these tests hold it to that BY REPRODUCING THE CONDITION: the pre-guard matcher below is argparse's own
earlier implementation, and it is installed on `argparse.ArgumentParser` for the duration of a test.
The control matters as much as the check -- with the fix removed and nothing else changed, the same
parser produces the CI message word for word, which is what makes the fix the thing that fixed it.
"""

from __future__ import annotations

import argparse
import re

import pytest

from genomeos import cli

#: What the CI message was, word for word (ci run at 2355bad, tests/test_roadmap_work.py:174).
CI_MESSAGE = "unrecognized arguments: a task"

#: The two orders. The first is what `cmd_work`'s own error message tells a lane to use and what the
#: failing test passes; the second is what the lane briefs use. Both must work, on every interpreter.
ORDERS = (
    ["work", "start", "--who", "lane-in-a-worktree", "a task"],
    ["work", "start", "a task", "--who", "lane-in-a-worktree"],
)


def pre_guard_match_arguments_partial(self, actions, arg_strings_pattern):  # noqa: N805 - a method
    """argparse's positional matcher BEFORE the guard, which is the behaviour that reddened CI.

    Measured against this parser shape on 2026-10-03: Python 3.9.6 gives `action='start'`, `task=None`
    and leaves `'a task'` over; Python 3.12.13 gives `task='a task'`. The two differ in exactly the
    three lines this function lacks, so installing it here reproduces the older interpreter without
    needing one installed.
    """
    result = []
    for i in range(len(actions), 0, -1):
        pattern = "".join(self._get_nargs_pattern(action) for action in actions[:i])
        match = re.match(pattern, arg_strings_pattern)
        if match is not None:
            result.extend([len(string) for string in match.groups()])
            break
    return result


@pytest.fixture
def pre_guard(monkeypatch):
    """The interpreter CI had: argparse's matcher without the guard, for the whole parser tree."""
    monkeypatch.setattr(
        argparse.ArgumentParser, "_match_arguments_partial", pre_guard_match_arguments_partial
    )


def test_argparse_still_has_the_private_names_the_fix_is_written_against():
    """If CPython renames either, the override stops being called and this says so instead of drifting.

    The failure mode is benign -- an interpreter that renamed them is long past the guard, so the parse
    would be correct anyway -- but a shim nothing checks is a shim nobody knows is dead.
    """
    for name in ("_match_arguments_partial", "_get_nargs_pattern"):
        assert hasattr(argparse.ArgumentParser, name), name
    assert "_match_arguments_partial" in vars(cli.StableParser), "the override is on the class"


def test_every_subcommand_s_parser_is_the_project_s_own_class():
    """`add_subparsers` defaults `parser_class` to `type(self)`; the fix is useless if that changes.

    The defect was in a SUBPARSER (`work`), not in the top-level parser, so this is the load-bearing
    half: checked against `work` by name and against every other subcommand too.
    """
    ap = cli.build_parser()
    assert isinstance(ap, cli.StableParser)
    actions = [a for a in ap._actions if isinstance(a, argparse._SubParsersAction)]
    assert len(actions) == 1, "one subparsers group"
    assert "work" in actions[0].choices, "the subcommand the defect was found in"
    assert all(isinstance(p, cli.StableParser) for p in actions[0].choices.values())


@pytest.mark.parametrize("argv", ORDERS, ids=["flags-before-the-task", "task-before-the-flags"])
def test_both_orders_reach_the_task_on_this_interpreter(argv):
    """Today's CPython, unmodified: the figure under test is `task`, and it is the task string."""
    args = cli.build_parser().parse_args(argv)
    assert (args.cmd, args.action, args.task, args.who) == (
        "work",
        "start",
        "a task",
        "lane-in-a-worktree",
    )


@pytest.mark.parametrize("argv", ORDERS, ids=["flags-before-the-task", "task-before-the-flags"])
def test_both_orders_reach_the_task_under_the_argparse_that_reddened_ci(argv, pre_guard):
    """The same two parses under the pre-guard matcher: the condition that broke CI, and they hold."""
    args = cli.build_parser().parse_args(argv)
    assert (args.cmd, args.action, args.task, args.who) == (
        "work",
        "start",
        "a task",
        "lane-in-a-worktree",
    )


def test_the_control_without_the_fix_reproduces_the_ci_message_word_for_word(monkeypatch, capsys, pre_guard):
    """THE CONTROL. The same parser with only the override removed, under the same matcher, fails.

    Without this the two tests above prove nothing: a parse that would have worked anyway cannot show
    what the fix is for. `delattr` on the subclass lets the pre-guard matcher through to every parser in
    the tree, including `work`'s, and nothing else about the parser changes.
    """
    monkeypatch.delattr(cli.StableParser, "_match_arguments_partial")
    with pytest.raises(SystemExit) as exit_info:
        cli.build_parser().parse_args(ORDERS[0])
    assert exit_info.value.code == 2
    assert CI_MESSAGE in capsys.readouterr().err
    # and the order the lane briefs use was never the broken one, which is why this was missed
    args = cli.build_parser().parse_args(ORDERS[1])
    assert args.task == "a task"
