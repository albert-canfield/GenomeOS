"""The roadmap as data: docs/ROADMAP.md read into areas, planned steps and finished items.

The Progress tab shows what is planned and what is done without anyone restating it, so the
roadmap stays the single plan and this module only reads it. An area of section 3 is a list of
bold-led bullets: the fixed keys (Goal, Code, Design, Data, Requirements, Missing, Next, Owner)
describe it, and every other bullet is a finished piece of work, usually dated in its lead.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

META = ("Goal", "Code", "Design", "Data", "Requirements", "Missing", "Next", "Owner")
_AREA = re.compile(r"^### ([A-Z])\. (.+)$")
_BULLET = re.compile(r"^- \*\*(.+?)\*\*(.*)$")
_DATE = re.compile(r"\b(20\d\d-\d\d-\d\d)\b")
# a numbered step inside a Next paragraph: "1. ", not a decimal ("0.4") nor a version ("v0.3")
_STEP = re.compile(r"(?:^|(?<=\s))(\d{1,2})\.\s+(?=\S)")


def _sections(text: str) -> dict[str, list[str]]:
    """Top-level sections (## N. Title) with their lines."""
    out: dict[str, list[str]] = {}
    cur = None
    for line in text.splitlines():
        m = re.match(r"^## (\d+)\. ", line)
        if m:
            cur = m.group(1)
            out[cur] = []
        elif cur:
            out[cur].append(line)
    return out


def _bullets(lines: list[str]) -> list[tuple[str, str]]:
    """Bold-led bullets as (lead, body), continuation lines and nested tables folded into the body."""
    items: list[tuple[str, list[str]]] = []
    for line in lines:
        m = _BULLET.match(line)
        if m:
            items.append((m.group(1).strip(), [m.group(2).strip()]))
        elif items and (line.startswith("  ") or not line.strip()):
            items[-1][1].append(line.strip())
    out = []
    for lead, raw in items:
        # paragraphs are separated by blank lines, table rows stay one per line, prose lines are joined
        body, prev = "", ""
        for s in raw:
            if not s:
                body += "\n\n" if body and not body.endswith("\n\n") else ""
            elif s.startswith("|"):
                sep = "\n" if prev.startswith("|") else ("\n\n" if body and not body.endswith("\n\n") else "")
                body += sep + s
            else:
                sep = "\n\n" if prev.startswith("|") else (" " if body and not body.endswith("\n\n") else "")
                body += sep + s
            prev = s or prev
        out.append((lead, body.strip()))
    return out


def _first_sentence(body: str, limit: int = 240) -> str:
    flat = re.sub(r"\s+", " ", body.split("\n|")[0]).strip()
    m = re.search(r"(?<=[a-z0-9)`%])\.\s+(?=[A-Z`*(])", flat)
    s = flat[: m.start() + 1] if m else flat
    return s if len(s) <= limit else s[: limit - 1].rstrip() + "…"


def split_steps(body: str) -> list[str]:
    """A Next paragraph ("1. a. 2. b.") as its steps, in order; a paragraph without numbers is one step."""
    flat = re.sub(r"\s+", " ", body).strip()
    marks = list(_STEP.finditer(flat))
    # keep only the running sequence 1, 2, 3 ... so a count inside a step is not a step
    seq, want = [], 1
    for m in marks:
        if int(m.group(1)) == want:
            seq.append(m)
            want += 1
    if not seq:
        return [flat] if flat else []
    return [
        flat[m.end() : (seq[i + 1].start() if i + 1 < len(seq) else len(flat))].strip()
        for i, m in enumerate(seq)
    ]


def step_state(step: str) -> str:
    """done | partial | blocked | planned, from how the roadmap words the step."""
    low = step.lower()
    if low.startswith("done") or low.startswith("closed"):
        return "partial" if re.search(r"\bnext\b|still to|remains", low) else "done"
    if "blocked" in low or "waits for" in low or "waits on" in low:
        return "blocked"
    if re.search(r"\bbuilt\b|\brunning\b|\blanded\b", low):
        return "partial"
    return "planned"


def _owner_names(body: str) -> list[str]:
    return sorted(set(re.findall(r"genomeos-[0-9a-z]{2}\b", body)))


def parse_areas(text: str) -> list[dict[str, Any]]:
    lines = _sections(text).get("3", [])
    areas: list[dict[str, Any]] = []
    chunk: list[str] = []
    head = None

    def close() -> None:
        if head is None:
            return
        letter, title = head
        area: dict[str, Any] = {"letter": letter, "title": title, "next": [], "done": [], "missing": ""}
        area["_text"] = "\n".join(chunk)
        for lead, body in _bullets(chunk):
            key = lead.rstrip(".").split(" ")[0].rstrip(",.")
            date = _DATE.search(lead)
            if key in META and not date or lead.rstrip(".") in META:
                if key == "Next":
                    area["next"] = [{"text": s, "state": step_state(s)} for s in split_steps(body)]
                elif key == "Missing":
                    area["missing"] = re.sub(r"\s+", " ", body).strip()
                elif key == "Owner":
                    area["owner"] = re.sub(r"\s+", " ", body).strip()
                    area["owners"] = _owner_names(body)
                elif key == "Code":
                    area["code"] = re.sub(r"\s+", " ", body).strip()
                else:
                    area[key.lower()] = re.sub(r"\s+", " ", body).strip()
            elif key == "Missing":  # "Missing, and blocked on data (checked 2026-09-12)." is a Missing note
                area["missing"] = (area["missing"] + " " + re.sub(r"\s+", " ", body)).strip()
            else:
                area["done"].append(
                    {
                        "title": re.sub(r"\s*\(\d{4}-\d\d-\d\d[^)]*\)\.?$", "", lead).rstrip(".").strip(),
                        "date": date.group(1) if date else None,
                        "summary": _first_sentence(body),
                        "body": body,
                    }
                )
        # "Missing. ... Done 2026-09-11: a, b and c." carries finished work inside the Missing paragraph
        m = re.search(r"\s*Done (20\d\d-\d\d-\d\d):\s*(.+)$", area["missing"])
        if m:
            area["missing"] = area["missing"][: m.start()].strip()
            area["done"].insert(
                0,
                {
                    "title": "Done since the review",
                    "date": m.group(1),
                    "summary": _first_sentence(m.group(2)),
                    "body": m.group(2),
                },
            )
        area["counts"] = {
            "done": len(area["done"]),
            "next": len(area["next"]),
            "blocked": sum(1 for s in area["next"] if s["state"] == "blocked"),
        }
        areas.append(area)

    for line in lines:
        m = _AREA.match(line)
        if m:
            close()
            head, chunk = (m.group(1), m.group(2).strip()), []
        elif head is not None:
            chunk.append(line)
    close()
    return areas


def _table(lines: list[str]) -> list[list[str]]:
    rows = []
    for line in lines:
        if not line.startswith("|") or set(line.replace("|", "").strip()) <= {"-", " "}:
            continue
        rows.append([c.strip() for c in line.strip().strip("|").split("|")])
    return rows[1:] if rows else []  # drop the header


def parse_milestones(text: str) -> list[dict[str, Any]]:
    out = []
    for cells in _table(_sections(text).get("6", [])):
        if len(cells) < 3:
            continue
        version = cells[0]
        state = (
            "done"
            if "✅" in version or cells[2].startswith("done")
            else "partial"
            if "◐" in version or "◑" in version
            else "planned"
        )
        name = re.sub(r"[*✅◐◑]", "", version).strip()
        out.append({"milestone": name, "goal": cells[1], "proof": cells[2], "state": state})
    return out


def parse_data_jobs(text: str) -> list[dict[str, Any]]:
    """Section 4's two tables: the per-chromosome layers and the one-off jobs, each with a done flag."""
    lines = _sections(text).get("4", [])
    out = []
    for cells in _table(lines):
        if cells[0] in ("Layer", "Job"):
            continue
        if len(cells) == 4:  # Layer | Done | Job | Notes
            done = (
                cells[2].strip() == "done"
                or bool(re.match(r"(\d+) of \1\b", cells[1]))
                or "closed" in cells[1]
            )
            out.append(
                {
                    "job": cells[0],
                    "state": "done" if done else "planned",
                    "detail": cells[1] + " — " + cells[3],
                }
            )
        elif len(cells) == 3:  # Job | Needs | Owner
            done = "**done" in cells[1] or cells[1].startswith("done")
            out.append(
                {
                    "job": cells[0],
                    "state": "done" if done else "planned",
                    "detail": cells[1].replace("**", ""),
                    "owner": cells[2],
                }
            )
    return out


def load(root: Path) -> dict[str, Any]:
    p = root / "docs" / "ROADMAP.md"
    if not p.exists():
        return {"areas": [], "milestones": [], "data_jobs": [], "reviewed": None}
    text = p.read_text()
    reviewed = _DATE.search(text[:600])
    areas = [{k: v for k, v in a.items() if k != "_text"} for a in parse_areas(text)]
    return {
        "reviewed": reviewed.group(1) if reviewed else None,
        "modified": p.stat().st_mtime,
        "areas": areas,
        "milestones": parse_milestones(text),
        "data_jobs": parse_data_jobs(text),
        "totals": {
            "done": sum(a["counts"]["done"] for a in areas),
            "next": sum(a["counts"]["next"] for a in areas),
            "blocked": sum(a["counts"]["blocked"] for a in areas),
        },
    }


# where a file lives says its area when the roadmap does not name it
_DIR_AREA = (
    ("genomeos/attribution/", "I"),
    ("genomeos/knowledge/homology", "J"),
    ("genomeos/knowledge/profiling", "J"),
    ("genomeos/knowledge/across", "J"),
    ("genomeos/decompile", "J"),
    ("genomeos/genome/", "B"),
    ("genomeos/predict/", "B"),
    ("genomeos/molecules/", "C"),
    ("genomeos/flow/", "C"),
    ("genomeos/twin/", "D"),
    ("data/twins/", "D"),
    ("genomeos/organism/", "E"),
    ("data/organisms/", "E"),
    ("genomeos/cancer/", "F"),
    ("genomeos/therapeutics/", "F"),
    ("genomeos/forge/", "H"),
    ("genomeos/design", "H"),
    ("genomeos/lang/", "A"),
    ("genomeos/ir/", "A"),
    ("genomeos/runtime/", "A"),
    ("genomeos/std/", "A"),
    ("genomeos/bio.py", "A"),
    ("genomeos/web/", "G"),
    ("genomeos/cli.py", "G"),
    ("genomeos/jobs.py", "G"),
)


def area_of_path(path: str, areas: list[dict[str, Any]]) -> str | None:
    """The area a changed file belongs to: the area whose text names it most, else where it lives."""
    p = Path(path)
    if p.parts[:1] == ("tests",) and p.stem.startswith("test_"):  # a test belongs where its module does
        name = f"/{p.stem[5:]}.py"
        counts = [(a.get("_text", "").count(name), a["letter"]) for a in areas]
        n, letter = max(counts, default=(0, None))
        return letter if n else None
    tail = f"{p.parent.name}/{p.name}" if p.parent.name else p.name
    brace = re.compile(rf"{re.escape(p.parent.name)}/\{{[^}}]*\b{re.escape(p.stem)}\b[^}}]*\}}")
    best, hits = None, 0
    for a in areas:
        text = a.get("_text") or a.get("code", "")
        n = text.count(tail) + len(brace.findall(text)) if p.parent.name else text.count(f"`{p.name}`")
        if n > hits:
            best, hits = a["letter"], n
    if best:
        return best
    for prefix, letter in _DIR_AREA:
        if path.startswith(prefix):
            return letter
    return None
