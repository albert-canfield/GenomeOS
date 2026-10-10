# SPDX-License-Identifier: AGPL-3.0-or-later
"""The worm embryo's MEASURED cell-cell contacts, distilled for the Body's `contacts:` table (area E).

data/organisms/celegans/contacts_embryo.tsv is curated: its own header says no measured contact table
was held locally, so every row was written from the published cleavage-stage geometry, and the facts the
founder signals turn on (ABa does not touch P2; ABp and EMS do) were taken from the same papers that
described the signals. This script replaces that borrowing with a measurement made for another purpose.

THE SOURCE. Cao J, Guan G, Ho VWS, et al. (2020) "Establishment of a morphological atlas of the
Caenorhabditis elegans embryo using deep-learning-based 4D segmentation", Nat Commun 11:6254,
doi:10.1038/s41467-020-19863-x (CShaper). Membrane-marker confocal movies of 17 wild-type embryos,
segmented and lineage-traced from the 4- to the 350-cell stage. Supplementary Data 5 lists the cell-cell
contact pairs at four successive filters:

    Sheet1  3,011 pairs  contact area > 0 at one time point or more, in one embryo or more
    Sheet2  2,790 pairs  and contact area >= 1/48 of either cell's surface
    Sheet3  2,276 pairs  and lasting >= 2 consecutive time points (about 3 min)
    Sheet4    479 pairs  and reproducible in all 17 embryos ("effective contacts")

Sheet4 is restricted to the 222 cells segmented at every time point in all 17 embryos. Article and
figshare record (doi:10.6084/m9.figshare.12839315) are CC BY 4.0. The file is 461,824 bytes; its
sha256 is pinned below and a different download is refused.

WHAT THE DISTILLED TABLES SAY, AND WHAT THEY DO NOT. A pair is listed without a time: Supplementary
Data 5 states that two cells touched, not when or for how long. Every row is therefore written at time
0, and the Body's own aliveness filter (`Body.neighbours`) does the rest: a cell's neighbours at t are
its measured partners that are alive at t. So a pair is treated as touching for the whole of the two
cells' common lifetime, which overstates a contact that held for part of it. Every area is 1.0: the
table is measured topology, not measured area (Supplementary Data 6 gives areas for ten known
signalling pairs only). The names are Sulston's, written with a trailing apostrophe in the
spreadsheet; the apostrophe is dropped and every one of the 286 names is checked against the reference
lineage, so a name the program does not grow cannot slip in unnoticed.

Two tables are written:

    contacts_embryo_cshaper.tsv           Sheet4, the 479 effective contacts (the program's table)
    contacts_embryo_cshaper_detected.tsv  Sheet1, the 3,011 detected pairs (the loosest reading,
                                          kept so a founder result can be checked against a contact
                                          seen once in one embryo)

    uv run --with xlrd python scripts/celegans_cshaper_contacts.py

xlrd reads the legacy .xls format and is not a project dependency, hence `--with`. No model calls. The
raw file is kept under data/knowledge/celegans/cshaper/ (git-ignored); the two tables are committed.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos.organism.reference import ReferenceLineage  # noqa: E402

URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19863-x/MediaObjects/"
    "41467_2020_19863_MOESM6_ESM.xls"
)
RAW = ROOT / "data/knowledge/celegans/cshaper/41467_2020_19863_MOESM6_ESM.xls"
SHA256 = "519c34349af8e8f6cbd675c665f9a880c99fe1ad482b469067f4a6078a855c99"
BYTES = 461_824
OUT = {
    "Sheet4": ROOT / "data/organisms/celegans/contacts_embryo_cshaper.tsv",
    "Sheet1": ROOT / "data/organisms/celegans/contacts_embryo_cshaper_detected.tsv",
}
EXPECTED_PAIRS = {"Sheet4": 479, "Sheet1": 3011}
FILTER = {
    "Sheet4": "the 479 effective contacts: area >= 1/48 of either cell's surface, >= 2 consecutive time "
    "points (about 3 min), reproducible in all 17 wild-type embryos (Samples 04-20)",
    "Sheet1": "the 3,011 detected pairs: contact area > 0 at one time point or more in one embryo or more "
    "of the 17",
}

HEADER = """\
# MEASURED cell-cell contacts of the early C. elegans embryo, for BioLang v0.4 §7.4 (`contacts:`).
# Rows are `time cell cell area`. Distilled by scripts/celegans_cshaper_contacts.py; do not edit by hand.
#
# Source: Cao J, Guan G, Ho VWS, et al. (2020) Nat Commun 11:6254, doi:10.1038/s41467-020-19863-x
# (CShaper), Supplementary Data 5, {sheet}: {filter}.
# File: {url}
# {bytes} bytes, sha256 {sha}. Licence: CC BY 4.0 (article and figshare doi:10.6084/m9.figshare.12839315).
# Segmented and lineage-traced from the 4- to the 350-cell stage; {pairs} pairs over {cells} cells, every
# name checked against the reference lineage.
#
# Every row is at time 0 and every area is 1.0, because the source states THAT two cells touched and not
# WHEN or HOW MUCH: the Body's aliveness filter makes a pair touch for the whole common lifetime of the two
# cells, which overstates a contact that held for part of it. This is measured topology, not measured
# area. Nothing before the 4-cell stage is in the source (AB-P1 is not listed).
#
# time\tcell\tcell\tarea
"""


def fetch() -> bytes:
    if not RAW.exists():
        RAW.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL, timeout=120) as r:  # noqa: S310 -- a fixed https URL
            RAW.write_bytes(r.read())
    data = RAW.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SHA256 or len(data) != BYTES:
        raise SystemExit(f"{RAW}: sha256 {digest} ({len(data)} bytes), pinned {SHA256} ({BYTES}): refused")
    return data


def pairs(data: bytes, sheet: str) -> list[tuple[str, str]]:
    import xlrd  # noqa: PLC0415 -- only this distil step needs it

    s = xlrd.open_workbook(file_contents=data).sheet_by_name(sheet)
    out: set[tuple[str, str]] = set()
    for r in range(2, s.nrows):  # row 0 is the title, row 1 the column names
        a, b = (str(v).strip().rstrip("'") for v in s.row_values(r)[:2])
        if a and b:
            out.add(tuple(sorted((a, b))))
    return sorted(out)


def main() -> None:
    data = fetch()
    ref = ReferenceLineage.load()
    for sheet, path in OUT.items():
        rows = pairs(data, sheet)
        cells = {c for p in rows for c in p}
        unknown = sorted(c for c in cells if c not in ref.cells)
        if unknown or len(rows) != EXPECTED_PAIRS[sheet]:
            raise SystemExit(f"{sheet}: {len(rows)} pairs, names not in the reference: {unknown[:10]}")
        head = HEADER.format(
            sheet=sheet,
            filter=FILTER[sheet],
            url=URL,
            bytes=f"{BYTES:,}",
            sha=SHA256,
            pairs=len(rows),
            cells=len(cells),
        )
        path.write_text(head + "".join(f"0\t{a}\t{b}\t1.0\n" for a, b in rows))
        print(f"{path.relative_to(ROOT)}: {len(rows)} pairs over {len(cells)} cells")


if __name__ == "__main__":
    main()
