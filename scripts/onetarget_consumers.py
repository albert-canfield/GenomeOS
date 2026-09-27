# SPDX-License-Identifier: AGPL-3.0-or-later
"""The last one-target consumers, read from the sweep's per-element cache beside the compact table.

    uv run python scripts/onetarget_consumers.py [--no-save]

Four modules read one target gene per element from the compact `all_elements` tables: the motif
transfer's sweep lead (`motif_transfer.tier_blocks`), the syntax tiling's `moved`, the syntax
candidates' node context (`candidates.node_context`) and the locus benchmark's gene-input layer
(`loci.read_gene_input`). Each is read here twice -- as it reads today and through
`attribution.targets.ElementResponses`, where every gene in the scorer's window is kept -- and the
registration below, committed before any of this was run, says which figures must come back to the
digit, which may move and in which direction, and what would falsify it. 0 AlphaGenome requests: the
cache is on disk, and a chromosome's archive is held one at a time.
"""

from __future__ import annotations

from genomeos.predict.enhancer_target import MIN_EFFECT

REGISTERED = "2026-09-27"

# A gene is "at the bar" at an element when the sweep's larger predicted change for it, the fall or the
# rise over all 371 tracks, reaches MIN_EFFECT. That is `predict_target`'s own rule applied to every
# gene of the window instead of to its head only, so an element has a gene at the bar exactly when the
# compact table names one; what the window adds is the other genes, never a different yes or no.
BAR = MIN_EFFECT

PRE_REGISTRATION = {
    "registered": REGISTERED,
    "one_registration_for": [
        "motif_transfer.tier_blocks (the sweep lead per block)",
        "syntax_tiling.moved (the tiling's own scored windows)",
        "candidates.node_context (node_targets_named, deleted_elements_on_block)",
        "loci.read_gene_input (the gene-input layer of the locus benchmark)",
    ],
    "shared_premise": (
        "The roadmap line says an element that moves a non-top gene reads as 'did not move'. It is "
        "registered here as expected FALSE: the compact table's head is the maximum over the window, so "
        "'some gene in the window is at the bar' and 'the head is at the bar' are the same statement. "
        "Measured on chr21 during the survey: 7,846 elements yes/yes, 4,293 no/no, 0 disagreements. "
        "What a one-target reader loses is the count and identity of the other genes, not the yes/no."
    ),
    "shared_controls": [
        "every element read from both sources: the window's head at the bar names the compact table's "
        "`predicted` gene and its coding head names `predicted_coding` -- 0 disagreements required",
        "an element the cache does not hold keeps its compact reading and is counted by name "
        "(not_cached), never scored as zero",
    ],
    "motif_transfer": {
        "reproduce_exactly": {
            "real_unknown.with_a_tested_element": 531,
            "real_unknown.with_a_sweep_lead": 331,
            "neutral.with_a_tested_element": 1181,
            "neutral.with_a_sweep_lead": 768,
            "by_case.with_a_sweep_lead": {
                "recent": 11,
                "relaxed": 167,
                "syntax": 19,
                "tolerant": 131,
                "unmeasured": 3,
            },
        },
        "reproduce_by_construction": (
            "the agreement AUROC 0.5422 [0.4937, 0.5934]: its label is with_a_sweep_lead per block, so if "
            "every block's label reproduces the AUROC cannot move; the count model is not re-run"
        ),
        "may_change": (
            "new field only, genes at the bar per block over its elements; expected >= the distinct head "
            "genes of the block, about 2.3 genes per moving element (chr21 survey 2.34)"
        ),
        "falsifier": "any block whose lead or tested count differs from the stored result",
    },
    "syntax_tiling": {
        "reproduce_exactly": "nothing: the registered run was cancelled unrun and has no result",
        "not_rerun": "scoring the tiles is model requests, and the key is held elsewhere",
        "change": (
            "`moved` is unchanged; the tiling's rows gain the genes at the bar read from the cache the "
            "scorer writes, so its result will not be one-target when it runs"
        ),
        "falsifier": "a cached element where `moved` and 'some gene at the bar' disagree (unit test)",
    },
    "candidates": {
        "reproduce_exactly": (
            "`node_targets_named`, `deleted_elements_in_node` and `deleted_elements_on_block` of the one "
            "candidate on chr22, the only candidate chromosome whose whole-chromosome run was complete "
            "when syntax_candidates_genome_wide was written on 2026-09-13"
        ),
        "expected_to_differ_from_the_stored_file_by_input": (
            "the other 68, whose chromosomes completed after 2026-09-13: the old reading recomputed today "
            "is the control for the window reading, and its difference from the stored file is attributed "
            "to the input and counted, not to the reader (the confound lane-scoring measured on 7f7c8c1)"
        ),
        "may_change": {
            "node_targets_named_window": (
                "for every gene, window count >= table count (a gene keeps every vote it had as head and "
                "gains one wherever it is at the bar below the head); strictly monotone, any gene whose "
                "window count is lower falsifies the reader"
            ),
            "top_node_gene": "two-sided, no direction registered: counted how many of the 69 change",
        },
        "not_rerun": "the reading, labels and set tests of syntax_candidates_genome_wide stand as stored",
    },
    "loci_gene_input": {
        "reproduce_exactly": (
            "every stored gene_input field of the 17 loci in loci_benchmark (2026-09-17, after the sweep "
            "finished): elements_scored, naming_a_coding_gene, the top-8 genes with their sums, target and "
            "rank_of_first_published_target; stored strict hits 10 of 17 (SHH_ZRS, HERC2_OCA2, MCM6_LCT, "
            "HBB_LCR, FTO_IRX3, BCL11A_enhancer, ABO, APP, HOXD, SOX9_PierreRobin), among 15 of 17"
        ),
        "window_rule": (
            "each element credits |effect| to every coding gene at the bar in its window, not only to its "
            "coding head; elements the cache does not hold (stated intervals, VISTA and lentiMPRA rows) "
            "keep their compact credit"
        ),
        "expected_direction": (
            "strict hits hold or rise (>= 10 of 17): a published target that is second at many elements "
            "gains credit (MYC_8q24 rank 2, PMP22_CMT1A rank 2 are the likeliest to flip); the median "
            "rank of the first published target over the 7 misses holds or improves"
        ),
        "falsifier": (
            "strict hits fall below 10, or any of the four loci gene_input carries alone in the headline "
            "(SHH_ZRS, HERC2_OCA2, FTO_IRX3, SOX9_PierreRobin) loses its hit: crediting every gene then "
            "crowds the published target, and the one-target sum was the better reader"
        ),
        "headline": (
            "target_derived 15/17 is NOT recomputed in place: swapping a layer of the shared benchmark is "
            "the owner's call. The counterfactual rate with the window layer in place of gene_input is "
            "reported beside it"
        ),
    },
}
