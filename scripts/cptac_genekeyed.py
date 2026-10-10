# SPDX-License-Identifier: AGPL-3.0-or-later
"""The seven gene-keyed CPTAC phospho profiles: gene keys and positions resolved by a registered rule,
or retired.

    uv run python scripts/cptac_genekeyed.py inputs     # HGNC set present; GENCODE v50 translations built
    uv run python scripts/cptac_genekeyed.py register   # data/results/cptac_genekeyed_registration.json
    uv run python scripts/cptac_genekeyed.py apply      # data/results/cptac_genekeyed_result.json

The rule, the gates and the wording are ``genomeos.molecules.cptac_keys.REGISTRATION``, written and
committed before ``apply`` ran on any site. ``apply`` reads only cached inputs (the cBioPortal entity
lists and RefSeq-keyed mapping of 2026-09-28, the HGNC complete set, the GENCODE translations, the compiled
UniProt definitions) until a profile has passed every mapping gate; only then does it fetch the values it
lacks for that profile from cBioPortal, summarised in memory as on 2026-09-28 (tumours with a value, median
log2 ratio), never per tumour.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import cptac_phospho as cp  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.molecules import cptac_keys as ck  # noqa: E402
from genomeos.molecules import ptm  # noqa: E402

CACHE = cp.CACHE
HGNC_FILE = CACHE / "hgnc_complete_set.txt"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
TRANSLATIONS = CACHE / "gencode_v50_translations.json.gz"
REFERENCE = Path("data/reference")
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
INCOMPLETE_TAGS = {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"}
ENTRY = "scripts/cptac_genekeyed.py"
OWN_CODE = (
    "genomeos/molecules/cptac_keys.py",
    "scripts/cptac_genekeyed.py",
    "scripts/cptac_phospho.py",
    "genomeos/molecules/ptm.py",
    "tests/test_cptac_genekeyed.py",
)
VALUE_DOWNLOAD_CAP = 300 * 1024**2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---- inputs -------------------------------------------------------------------------------------------


def _chrom_sequence(chrom: str) -> str:
    parts = []
    with gzip.open(REFERENCE / f"{chrom}.fa.gz", "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def build_translations() -> dict[str, Any]:
    """ENSG -> distinct GENCODE v50 translations with a complete CDS, and every exclusion counted."""
    from genomeos.runtime.central_dogma import (
        SELENOCYSTEINE_CODE,
        STANDARD_CODE,
        VERTEBRATE_MITOCHONDRIAL_CODE,
        translate_cds,
    )

    out: dict[str, list[str]] = defaultdict(list)
    why: Counter = Counter()
    gff_digests = {}
    for chrom in CHROMS:
        gff = REFERENCE / f"gencode_v50_{chrom}.gff3.gz"
        gff_digests[gff.name] = sha256(gff)
        txs: dict[str, dict[str, Any]] = {}
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS" or f[0] != chrom:
                    continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                t = txs.setdefault(
                    a["transcript_id"],
                    {
                        "gene": a["gene_id"].split(".")[0],
                        "strand": f[6],
                        "tags": set(a.get("tag", "").split(",")) - {""},
                        "segs": [],
                    },
                )
                t["segs"].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] != "." else 0))
        seq = _chrom_sequence(chrom)
        table = VERTEBRATE_MITOCHONDRIAL_CODE if chrom == "chrM" else None
        for _tid, t in txs.items():
            if t["tags"] & INCOMPLETE_TAGS:
                why["incomplete CDS tag"] += 1
                continue
            dna = "".join(seq[s:e] for s, e, _ in sorted(t["segs"]))
            if t["strand"] == "-":
                dna = dna.translate(_COMPLEMENT)[::-1]
            first = max(t["segs"]) if t["strand"] == "-" else min(t["segs"])  # 5'-most segment
            dna = dna[first[2] :]
            if len(dna) % 3:
                why["CDS length not a multiple of three"] += 1
                continue
            tab = table or (SELENOCYSTEINE_CODE if "seleno" in t["tags"] else STANDARD_CODE)
            prot = translate_cds(dna.replace("T", "U"), tab)
            if len(prot) != len(dna) // 3 - 1:
                why["stop before the final codon, or no final stop"] += 1
                continue
            why["kept"] += 1
            if prot not in out[t["gene"]]:
                out[t["gene"]].append(prot)
        print(chrom, len(txs), dict(why), flush=True)
    meta = {
        "source": "GENCODE v50 per-chromosome GFF3 (data/reference/gencode_v50_chr*.gff3.gz) and GRCh38 "
        "chromosome FASTA (data/reference/chr*.fa.gz), CDS rows only",
        "gff3_sha256": gff_digests,
        "transcripts": dict(why),
        "genes": len(out),
        "distinct_translations": sum(len(v) for v in out.values()),
    }
    with gzip.open(TRANSLATIONS, "wt") as fh:
        json.dump({"meta": meta, "translations": out}, fh)
    return meta


def inputs() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    if not HGNC_FILE.exists():
        HGNC_FILE.write_bytes(cp._http(HGNC_URL))
    print("HGNC", HGNC_FILE, sha256(HGNC_FILE))
    if not TRANSLATIONS.exists():
        print(json.dumps(build_translations(), indent=1))
    else:
        print("translations present", TRANSLATIONS)


# ---- shared loaders -----------------------------------------------------------------------------------


def load_compiled() -> dict[str, dict[str, Any]]:
    """accession -> {gene, sequence, sites, ensg} from the compiled UniProt definitions."""
    out: dict[str, dict[str, Any]] = {}
    for p in sorted(ptm.CACHE.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" not in d:
            continue
        acc = ptm.accession_of(d)
        if not acc:
            continue
        ident = d["sections"].get("identity", {}).get("items", {})
        origin = (d["sections"].get("genomic_origin") or {}).get("items") or {}
        e = out.setdefault(
            acc,
            {
                "gene": ident.get("gene") or d.get("gene"),
                "sequence": ptm.sequence_of(d),
                "sites": [(s["start"], s["class"]) for s in ptm.sites(d) if isinstance(s.get("start"), int)],
                "ensg": None,
            },
        )
        if origin.get("gene_id") and not e["ensg"]:
            e["ensg"] = str(origin["gene_id"]).split(".")[0]
    return out


def load_translations() -> tuple[dict[str, list[str]], dict[str, Any]]:
    with gzip.open(TRANSLATIONS, "rt") as fh:
        d = json.load(fh)
    return d["translations"], d["meta"]


def input_digests() -> dict[str, Any]:
    metas = sorted(CACHE.glob("meta_*.json"))
    return {
        "hgnc_complete_set": {"path": str(HGNC_FILE), "sha256": sha256(HGNC_FILE), "url": HGNC_URL},
        "gencode_v50_translations": {"path": str(TRANSLATIONS), "sha256": sha256(TRANSLATIONS)},
        "refseq_keyed_mapping": {"path": str(CACHE / "mapped.json"), "sha256": sha256(CACHE / "mapped.json")},
        "cbioportal_entity_lists": {p.name: sha256(p) for p in metas},
    }


def manifest_for(extra_inputs: list[dict[str, Any]], parameters: dict[str, Any]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": f"cBioPortal:{p}",
                "version": "generic-assay entity lists, public API, fetched 2026-09-27/28",
                "licence": "ODbL-1.0",
            }
            for p in (*ck.GENE_KEYED_PROFILES, *ck.REFSEQ_KEYED_PROFILES)
        ]
        + [
            {"accession": "HGNC complete set", "version": "fetched 2026-10-10", "licence": "CC0 1.0"},
            {"accession": "GENCODE v50 (GRCh38.p14)", "version": "per-chromosome GFF3", "licence": "public"},
            {
                "accession": "UniProtKB/Swiss-Prot compiled definitions",
                "version": "2026_03",
                "licence": "CC BY 4.0",
            },
        ],
        "inputs": [
            mf.input_entry(str(CACHE), None),
            mf.input_entry(str(ptm.CACHE), None),
            *extra_inputs,
        ],
        "assembly": "n/a for the result: protein residue positions on UniProt canonical sequences; the "
        "candidate isoforms were translated from GRCh38 with GENCODE v50",
        "coordinates": {"base": 1, "interval": "closed", "note": "residue positions"},
        "code": {},
        "parameters": parameters,
        "exclusions": [
            "multi-site and unlocalised entities (parse stage of scripts/cptac_phospho.py)",
            "gene keys that are ambiguous or unresolved under the registered HGNC rule",
            "sites that are ambiguous, carried only off the canonical, or unresolved",
            "several entities on one canonical site in one profile",
            "profiles failing any registered gate",
        ],
        "partitions": {
            "calibration": "the RefSeq-keyed profiles' mapped sites, whose canonical position is known",
            "application": "the seven gene-keyed profiles",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


# ---- register -----------------------------------------------------------------------------------------


def register() -> None:
    from genomeos.results import save_result

    mf.trace_begin()
    tr, tmeta = load_translations()
    compiled = load_compiled()
    with_tr = [d for d in compiled.values() if d["ensg"] and d["ensg"] in tr]
    canon_in = sum(1 for d in with_tr if d["sequence"] in tr[d["ensg"]])
    payload = {
        "registered_by": "lane-cptac",
        "registered_on": "2026-10-10",
        "registration": {k: (list(v) if isinstance(v, tuple) else v) for k, v in ck.REGISTRATION.items()},
        "profiles": list(ck.GENE_KEYED_PROFILES),
        "calibration_profiles": list(ck.REFSEQ_KEYED_PROFILES),
        "reproduced_before_registration": {
            "gene_keyed_residue_mismatch": "31,424 of 219,035 single-site entities (0.1435)",
            "per_profile": {
                "brca_cptac_2020_phosphoproteome": "3,449 of 25,375 (0.1359)",
                "ucec_cptac_2020_phosphoproteome": "10,354 of 71,125 (0.1456)",
                "coad_cptac_2019_phosphoprotein_quantification": "1,292 of 30,331 (0.0426)",
                "brca_tcga_phosphoprotein_quantification": "9,922 of 56,244 (0.1764)",
                "ov_tcga_phosphoprotein_quantification": "3,094 of 16,075 (0.1925)",
                "brca_tcga_pan_can_atlas_2018_phosphoprotein_quantification": "2,594 of 15,993 (0.1622)",
                "ov_tcga_pan_can_atlas_2018_phosphoprotein_quantification": "719 of 3,892 (0.1847)",
            },
            "gene_key_drops": "8,085 single-site entities, 498 distinct symbols, none naming several "
            "compiled entries (all named none)",
            "note": "docs/PROTEIN.md gives 219,179: the same count before the 144 duplicate-site drops; "
            "both round to 14.3%",
        },
        "inputs_fixed_at_registration": input_digests(),
        "gencode_translations": tmeta,
        "input_checks": {
            "compiled_entries": len(compiled),
            "with_gencode_translations": len(with_tr),
            "uniprot_canonical_identical_to_a_gencode_translation": canon_in,
            "note": "an input check of the translation build, not an outcome of the rule",
        },
        "outcome_seen": False,
    }
    man = manifest_for(
        [mf.input_entry(str(HGNC_FILE), None), mf.input_entry(str(TRANSLATIONS), None)],
        {"flank": 7, "decoy_offsets": list(ck.REGISTRATION["decoy_offsets"])},
    )
    path = save_result("cptac_genekeyed_registration", payload, manifest=man)
    print(path)


# ---- apply --------------------------------------------------------------------------------------------


class Resolver:
    def __init__(self, compiled: dict[str, dict[str, Any]], hgnc: ck.Hgnc, tr: dict[str, list[str]]):
        self.compiled = compiled
        self.hgnc = hgnc
        self.tr = tr
        self.by_ens: dict[str, list[str]] = defaultdict(list)
        for a, d in compiled.items():
            if d["ensg"]:
                self.by_ens[d["ensg"]].append(a)
        self._keys: dict[str, tuple[str, str | None, str]] = {}

    def key(self, symbol: str) -> tuple[str, str | None, str]:
        if symbol not in self._keys:
            self._keys[symbol] = ck.resolve_gene_key(symbol, self.hgnc, set(self.compiled), self.by_ens)
        return self._keys[symbol]

    def candidates(self, symbol: str, acc: str) -> list[str]:
        ensg = self.compiled[acc]["ensg"]
        if not ensg:
            _, hid, _ = self.hgnc.record(symbol)
            ensg = self.hgnc.ensembl.get(hid or "")
        return self.tr.get(ensg or "", [])

    def site(self, symbol: str, acc: str, residue: str, position: int) -> tuple[str, int | None, str]:
        return ck.resolve_site(
            residue, position, self.compiled[acc]["sequence"], self.candidates(symbol, acc)
        )


def single_sites(profile: str) -> tuple[list[tuple[str, dict]], Counter]:
    rows, drops = [], Counter()
    for e in json.loads((CACHE / f"meta_{profile}.json").read_text()):
        r = ptm.parse_cptac_entity(profile, e)
        if "drop" in r:
            drops["parse: " + r["drop"]] += 1
        else:
            rows.append((e["stableId"], r))
    return rows, drops


def run_profile(profile: str, res: Resolver, old_by_sid: dict[str, tuple[str, int]]) -> dict[str, Any]:
    rows, drops = single_sites(profile)
    keys: dict[str, list[Any]] = {}
    key_entities: Counter = Counter()
    status: Counter = Counter()
    hows: Counter = Counter()
    out: list[tuple[str, str, int, str]] = []
    unresolved_sites: list[str] = []
    ambiguous_sites: list[str] = []
    decoy: Counter = Counter()
    vs_old: Counter = Counter()
    for sid, r in rows:
        k_status, acc, k_how = res.key(r["gene"])
        keys[r["gene"]] = [k_status, acc, k_how]
        key_entities[r["gene"]] += 1
        if acc is None:
            status["gene key " + k_status] += 1
            continue
        s, cpos, how = res.site(r["gene"], acc, r["residue"], r["position"])
        status[s] += 1
        if s == "resolved":
            hows[how] += 1
            out.append((sid, acc, cpos, r["residue"]))
        elif s == "unresolved":
            unresolved_sites.append(sid)
        elif s == "ambiguous":
            ambiguous_sites.append(sid)
        for off in ck.REGISTRATION["decoy_offsets"]:
            if r["position"] + off >= 1:
                decoy[res.site(r["gene"], acc, r["residue"], r["position"] + off)[0]] += 1
        if sid in old_by_sid:
            o = old_by_sid[sid]
            vs_old[
                "same canonical site"
                if s == "resolved" and (acc, cpos) == o
                else "another canonical site"
                if s == "resolved"
                else s
            ] += 1
    n = Counter((a, p) for _, a, p, _ in out)
    keep = [x for x in out if n[(x[1], x[2])] == 1]
    dup = len(out) - len(keep)
    with_key = sum(v for k, v in status.items() if not k.startswith("gene key"))
    u = status["unresolved"]
    w = ck.estimated_wrong(u, decoy["resolved"], decoy["unresolved"])
    resolved = status["resolved"]
    k_by = Counter(v[0] for v in keys.values())
    return {
        "entities": len(rows) + sum(drops.values()),
        "parse_drops": dict(drops.most_common()),
        "single_site_entities": len(rows),
        "gene_keys": {
            "distinct_symbols": len(keys),
            "resolved": k_by["resolved"],
            "ambiguous": k_by["ambiguous"],
            "unresolved": k_by["unresolved"],
            "resolved_how": dict(Counter(v[2] for v in keys.values() if v[0] == "resolved").most_common()),
            "ambiguous_listed": {
                g: [v[2], key_entities[g]] for g, v in sorted(keys.items()) if v[0] == "ambiguous"
            },
            "unresolved_listed": {
                g: [v[2], key_entities[g]] for g, v in sorted(keys.items()) if v[0] == "unresolved"
            },
            "entities_with_ambiguous_key": status["gene key ambiguous"],
            "entities_with_unresolved_key": status["gene key unresolved"],
        },
        "sites_with_resolved_key": with_key,
        "sites": {
            "resolved": resolved,
            "ambiguous": status["ambiguous"],
            "off_canonical": status["off_canonical"],
            "unresolved": u,
            "resolved_how": dict(hows.most_common()),
            "duplicate_canonical_site_dropped": dup,
            "kept": len(keep),
        },
        "unresolved_rate": round(u / with_key, 4) if with_key else None,
        "decoys": dict(decoy.most_common()),
        "estimated_wrong": round(w, 1) if w is not None else None,
        "estimated_wrong_rate": round(w / resolved, 4) if w is not None and resolved else None,
        "against_the_canonical_only_mapping_of_2026_09_28": dict(vs_old.most_common()),
        "unresolved_entities": sorted(unresolved_sites),
        "ambiguous_entities": sorted(ambiguous_sites),
        "_kept": keep,
    }


def calibrate(res: Resolver, mapped: dict[str, list[list[Any]]]) -> dict[str, Any]:
    """The rule on RefSeq-keyed sites whose canonical position is known from their named RefSeq protein."""
    per: dict[str, Any] = {}
    tot: Counter = Counter()
    for prof in ck.REFSEQ_KEYED_PROFILES:
        rows, _ = single_sites(prof)
        by_sid = {sid: r for sid, r in rows}
        c: Counter = Counter()
        decoy: Counter = Counter()
        for sid, acc, pos, _res, tier in mapped.get(prof, []):
            if tier != "refseq" or sid not in by_sid:
                continue
            r = by_sid[sid]
            k_status, kacc, _ = res.key(r["gene"])
            if kacc is None:
                c["gene key " + k_status] += 1
                continue
            if kacc != acc:
                c["gene key names another entry"] += 1
            s, cpos, _ = res.site(r["gene"], kacc, r["residue"], r["position"])
            if s == "resolved":
                c["correct" if (kacc, cpos) == (acc, pos) else "wrong"] += 1
            else:
                c[s] += 1
            for off in ck.REGISTRATION["decoy_offsets"]:
                if r["position"] + off >= 1:
                    decoy[res.site(r["gene"], kacc, r["residue"], r["position"] + off)[0]] += 1
        w = ck.estimated_wrong(c["unresolved"], decoy["resolved"], decoy["unresolved"])
        per[prof] = {
            **dict(c.most_common()),
            "wrong_rate": round(c["wrong"] / (c["correct"] + c["wrong"]), 4)
            if c["correct"] + c["wrong"]
            else None,
            "decoys": dict(decoy.most_common()),
            "estimated_wrong": round(w, 1) if w is not None else None,
        }
        tot.update(c)
    resolved = tot["correct"] + tot["wrong"]
    rate = round(tot["wrong"] / resolved, 4) if resolved else None
    bound = ck.REGISTRATION["gates"]["G0_rule_wrong_rate_on_refseq_truth_below"]
    return {
        "pooled": dict(tot.most_common()),
        "wrong_rate": rate,
        "G0_pass": rate is not None and rate < bound,
        "per_profile": per,
    }


def apply() -> None:
    from genomeos.results import load_result, save_result

    mf.trace_begin()
    reg = load_result("cptac_genekeyed_registration")
    if not reg or reg.get("registration", {}).get("gates") is None:
        raise SystemExit("no committed registration: run `register` and commit it first")
    gates = ck.REGISTRATION["gates"]
    compiled = load_compiled()
    hgnc = ck.Hgnc.from_tsv(HGNC_FILE)
    tr, tmeta = load_translations()
    res = Resolver(compiled, hgnc, tr)
    m = json.loads((CACHE / "mapped.json").read_text())
    mapped = m["mapped"]
    cal = calibrate(res, mapped)
    print("calibration", json.dumps({k: v for k, v in cal.items() if k != "per_profile"}), flush=True)

    curated_phospho = {(a, p) for a, d in compiled.items() for p, cls in d["sites"] if cls == "phospho"}
    lo, hi = gates["G1_median_of_site_medians_within"]
    client = None
    per_profile: dict[str, Any] = {}
    committed: list[str] = []
    joined: dict[tuple[str, int], dict[str, Any]] = {}
    wording: dict[str, str] = {}
    for prof in ck.GENE_KEYED_PROFILES:
        old = {sid: (acc, pos) for sid, acc, pos, _r, tier in mapped.get(prof, []) if tier == "gene"}
        st = run_profile(prof, res, old)
        keep = st.pop("_kept")
        failed = []
        if not cal["G0_pass"]:
            failed.append(
                f"G0 rule wrong rate on the RefSeq-keyed truth {cal['wrong_rate']} not below "
                f"{gates['G0_rule_wrong_rate_on_refseq_truth_below']}"
            )
        if st["unresolved_rate"] is None or st["unresolved_rate"] >= gates["G2_unresolved_rate_below"]:
            failed.append(
                f"G2 unresolved {st['sites']['unresolved']:,} of {st['sites_with_resolved_key']:,} "
                f"({st['unresolved_rate']}) not below {gates['G2_unresolved_rate_below']}"
            )
        if (
            st["estimated_wrong_rate"] is None
            or st["estimated_wrong_rate"] >= gates["G3_estimated_wrong_rate_below"]
        ):
            failed.append(
                f"G3 estimated wrong resolutions {st['estimated_wrong']} of {st['sites']['resolved']:,} "
                f"({st['estimated_wrong_rate']}) not below {gates['G3_estimated_wrong_rate_below']}"
            )
        hit = [x for x in keep if (x[1], x[2]) in curated_phospho]
        st["kept_on_curated_phospho_site"] = len(hit)
        vals: dict[str, tuple[int, float]] = {}
        if not failed:
            # values only for a profile that passed every mapping gate
            old_cache = CACHE / f"summary_{prof}.json"
            if old_cache.exists():
                vals.update({k: tuple(v) for k, v in json.loads(old_cache.read_text()).items()})
            new_cache = CACHE / f"summary_genekeyed_{prof}.json"
            if new_cache.exists():
                vals.update({k: tuple(v) for k, v in json.loads(new_cache.read_text()).items()})
            missing = [x[0] for x in hit if x[0] not in vals]
            st["values_fetched_for"] = len(missing)
            if missing:
                if client is None:
                    from genomeos.cancer.cbioportal import CBioPortal

                    client = CBioPortal()
                cp.DOWNLOAD_CAP = VALUE_DOWNLOAD_CAP
                got = cp.fetch_summaries(client, prof, missing)
                prior = json.loads(new_cache.read_text()) if new_cache.exists() else {}
                prior.update(got)
                new_cache.write_text(json.dumps(prior))
                vals.update(got)
            meds = [vals[x[0]][1] for x in hit if x[0] in vals]
            st["with_a_value"] = len(meds)
            pm = round(median(meds), 3) if meds else None
            st["median_of_site_medians"] = pm
            if pm is None or not lo <= pm <= hi:
                failed.append(f"G1 median of site medians {pm} outside [{lo}, {hi}]")
        st["gates_failed"] = failed
        if failed:
            st["outcome"] = "retired"
            wording[prof] = ck.REGISTRATION["wording_retired"].format(profile=prof, failed="; ".join(failed))
        else:
            st["outcome"] = "committed"
            committed.append(prof)
            wording[prof] = ck.REGISTRATION["wording_committed"].format(
                profile=prof,
                resolved=f"{st['sites']['resolved']:,}",
                sites=f"{st['sites_with_resolved_key']:,}",
                unresolved=f"{st['sites']['unresolved']:,}",
                unresolved_rate=f"{st['unresolved_rate']:.2%}",
                ambiguous=f"{st['sites']['ambiguous']:,}",
                estimated_wrong_rate=f"{st['estimated_wrong_rate']:.2%}",
            )
            for sid, acc, pos, residue in hit:
                if sid in vals:
                    j = joined.setdefault((acc, pos), {"residue": residue, "rows": []})
                    j["rows"].append(
                        [committed.index(prof), int(vals[sid][0]), round(float(vals[sid][1]), 3)]
                    )
        per_profile[prof] = st
        print(prof, st["outcome"], st["gates_failed"], flush=True)

    sites: dict[str, list[list[Any]]] = defaultdict(list)
    for (acc, pos), v in sorted(joined.items()):
        sites[acc].append([pos, v["residue"], sorted(v["rows"])])
    off = [
        f"{a}:{v['residue']}{p}" for (a, p), v in joined.items() if v["residue"] not in ptm.PHOSPHO_ACCEPTORS
    ]
    payload = {
        "field": ptm.CPTAC_FIELD,
        "meaning": ptm.CPTAC_MEANING,
        "tier": "gene-keyed: gene symbol resolved by the registered HGNC rule, position carried from the "
        "isoform that carries the residue onto the UniProt canonical (genomeos/molecules/cptac_keys.py)",
        "wording": wording,
        "committed_profiles": committed,
        "retired_profiles": [p for p in ck.GENE_KEYED_PROFILES if p not in committed],
        "calibration_on_refseq_keyed_truth": cal,
        "per_profile": per_profile,
        "joined_committed": len(joined),
        "proteins_with_joined_sites": len(sites),
        "curated_phospho_sites": len(curated_phospho),
        "off_acceptor_joined": off,
        "studies": committed,
        "sites_format": (
            "accession -> [[position, residue, [[study index, tumours with a value, median log2 ratio]]]]"
        ),
        "sites": dict(sites),
        "registration_result": "data/results/cptac_genekeyed_registration.json",
        "gencode_translations": tmeta,
        "licence": ptm.CPTAC_LICENCE.replace(
            "Residue mapping: UniProt", "Gene keys: HGNC (CC0). Residue mapping: GENCODE v50, UniProt"
        ),
        "downloaded_bytes": cp._downloaded,
    }
    man = manifest_for(
        [mf.input_entry("data/results/cptac_genekeyed_registration.json", None)],
        {
            "flank": 7,
            "decoy_offsets": list(ck.REGISTRATION["decoy_offsets"]),
            "gates": {k: (list(v) if isinstance(v, tuple) else v) for k, v in gates.items()},
            "value_download_cap_bytes": VALUE_DOWNLOAD_CAP,
        },
    )
    path = save_result("cptac_genekeyed_result", payload, manifest=man)
    print(json.dumps(wording, indent=1))
    print(path)


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    if stage == "inputs":
        inputs()
    elif stage == "register":
        register()
    elif stage == "apply":
        apply()
    else:
        raise SystemExit("usage: cptac_genekeyed.py inputs|register|apply")


if __name__ == "__main__":
    main()
