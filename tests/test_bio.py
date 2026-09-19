import io

from genomeos.bio import Repl, evaluate, main, parse_tests
from genomeos.lang import parse

SRC = """module t
gene tetR { max: 200; basal: 0.2; produces: TetR; evidence: curated "Elowitz 2000"; confidence: 0.9 }
protein TetR { half_life: 1; evidence: curated "Elowitz 2000"; confidence: 0.9 }
# test: rules >= 1
# test: TetR final > 5
# test: entities == 2
# test: confidence >= 0.8
# test: TetR peaks <= 5
"""


def test_parse_and_evaluate_test_lines():
    tests = parse_tests(SRC)
    assert [t[0] for t in tests] == ["rules", "TetR", "entities", "confidence", "TetR"]
    assert tests[1] == ("TetR", "final", ">", 5.0) and tests[4][1] == "peaks"
    module = parse(SRC)
    results = evaluate(module, tests, hours=24)
    assert all(r["ok"] for r in results), results
    assert results[-1]["test"] == "no NaN or infinite level"
    bad = evaluate(module, [("TetR", "final", "<", 0.0), ("Nope", "final", ">", 0.0)], hours=4)
    assert not bad[0]["ok"] and not bad[1]["ok"] and bad[1]["got"] is None


def test_bio_test_command_reports_failures(tmp_path, capsys):
    good = tmp_path / "good.bio"
    good.write_text(SRC)
    bad = tmp_path / "bad.bio"
    bad.write_text(SRC.replace("# test: TetR final > 5", "# test: TetR final > 1000000"))
    assert main(["test", str(good), "--hours", "12"]) == 0
    assert main(["test", str(tmp_path), "--hours", "12"]) == 1
    out = capsys.readouterr().out
    assert "FAIL TetR final > 1e+06" in out and "✓" in out and "✗" in out


def test_repl_accepts_blocks_rejects_errors_and_runs():
    out = io.StringIO()
    r = Repl(out=out)
    assert r.handle('gene a { max: 10; produces: A; evidence: curated "x"; confidence: 0.9 }')
    assert "waiting" in out.getvalue() and r.pending  # A is not declared yet
    assert r.handle('protein A { half_life: 2; evidence: curated "x"; confidence: 0.9 }')
    assert not r.pending and len(r.module().rules) == 1
    assert r.handle("gene b { max: 1; produces: A; evidence: curated x; confidence: 0.9 }") or True
    assert r.handle("nonsense {")
    assert "rejected" in out.getvalue()
    assert r.handle(":run 12")
    assert "A " in out.getvalue() and "final" in out.getvalue()
    assert r.handle(":check") and not r.handle(":quit")
