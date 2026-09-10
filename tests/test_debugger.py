"""Task 4.4: breakpoints on biological conditions and an evidence trace."""

from genomeos.lang import parse_file
from genomeos.runtime.debugger import AgeingDebugger, Breakpoint, NetworkDebugger


def test_breakpoint_parsing():
    bp = Breakpoint.parse("TetR > 50")
    assert (bp.variable, bp.op, bp.value) == ("TetR", ">", 50.0)
    assert Breakpoint.parse("senescent == true").value == 1.0


def test_network_breakpoint_and_explanation():
    dbg = NetworkDebugger(parse_file("data/demo/repressilator.bio"))
    dbg.set_initial({"TetR": 10.0})
    dbg.add_breakpoint("LacI > 30")
    hit = dbg.run_until_break(max_hours=60)
    assert hit is not None and dbg.session.state["LacI"] > 30
    lines = dbg.explain("lacI.mRNA")
    text = "\n".join(line.format() for line in lines)
    assert "CI inhibits lacI" in text and "Elowitz" in text  # the rule and its source appear
    assert any(line.confidence == 0.9 for line in lines)
    prot = dbg.explain("LacI")
    assert any("translation" in line.message for line in prot)


def test_ageing_breakpoint_and_event_trace():
    dbg = AgeingDebugger("hematopoietic_stem", seed=2)
    dbg.add_breakpoint("divisions >= 10")
    hit = dbg.step(years=80)
    assert hit is not None and dbg.session.state["divisions"] >= 10 and dbg.session.time < 80
    assert any(line.subject == "divide" and "Harley" in line.evidence for line in dbg.session.trace)
    lines = dbg.explain()
    assert any(line.subject == "telomere_loss_per_division_bp" and line.confidence == 0.6 for line in lines)
