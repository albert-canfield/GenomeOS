# SPDX-License-Identifier: AGPL-3.0-or-later
"""Text lines for `genomeos flow --tissue/--variant`: the tissue's isoform beside the canonical one.

Kept apart from genomeos/flow/tissue.py, whose bytes the registered census pins
(data/results/flow_tissue_census_registration.json)."""

from __future__ import annotations

import re
from typing import Any

from genomeos.flow.tissue import ESTABLISHED, TissueTrace

_VARIANT = re.compile(r"^(?:[\w.]+:)?([\d,]+)\s*([ACGTacgt])\s*>\s*([ACGTacgt])$")


def parse_variant(text: str) -> tuple[int, str, str]:
    """'25897620 G>A', '25897620G>A' or 'chr21:25,897,620 G>A' (1-based) → (0-based position, ref, alt)."""
    m = _VARIANT.match(text.strip())
    if not m:
        raise ValueError("variant format: POSITION REF>ALT, e.g. 25897620 G>A (1-based)")
    return int(m.group(1).replace(",", "")) - 1, m.group(2).upper(), m.group(3).upper()


def _call_text(call: dict[str, Any]) -> str:
    tpm = f"{call['tpm']} TPM" if call.get("tpm") is not None else "no TPM"
    nxt = f", next {call['next_transcript']} {call['next_tpm']} TPM" if call.get("next_transcript") else ""
    return f"{call.get('transcript') or '-'} {tpm}{nxt}"


def flow_lines(tt: TissueTrace, variant: tuple[int, str, str] | None = None) -> list[str]:
    """The lines `genomeos flow` prints for a tissue: the isoform and, with a variant, both consequences."""
    call = tt.call
    canon = f"canonical {tt.canonical.transcript} {len(tt.canonical.protein)} aa"
    if tt.status == ESTABLISHED:
        same = " (the canonical itself)" if tt.same_as_canonical else ""
        lines = [
            f"  tissue     {call['tissue']}: {tt.isoform.transcript} ({tt.isoform_id}){same}, "
            f"{len(tt.isoform.protein)} aa; {canon}",
            f"             dominant isoform {_call_text(call)}, share {call.get('share')}  "
            "[measured: GTEx v8 median, not this person's expression]",
        ]
    else:
        lines = [
            f"  tissue     {call['tissue']}: UNKNOWN ({call['reason']}: {_call_text(call)}); {canon}",
            "             the tissue's isoform is not established, so its consequence is UNKNOWN, "
            "never the canonical's",
        ]
    if variant is not None:
        pos, ref, alt = variant
        s = tt.substitute(pos, ref, alt)
        c = s["canonical"]
        lines.append(
            f"  variant    {pos + 1} {ref}>{alt}: in {call['tissue']} {s['consequence']} "
            f"{s.get('hgvs_p') or ''}".rstrip()
            + f"; on the canonical {c['consequence']} {c.get('hgvs_p') or ''}".rstrip()
            + (f" ({s['agreement']})" if s.get("agreement") else "")
            + "  [predicted by the local trace]"
        )
    return lines
