"""
mutation_scorer.py
-------------------
Core logic & bioinformatics engine for the Protein Sequence Mutation Scorer.

Owned by: Member 1 - Abhiram Sripathi (SE25UBIT014)
    - Sequence alignment (global / Needleman-Wunsch, via Biopython)
    - Mutation detection (substitutions, insertions, deletions)
    - Impact scoring (BLOSUM62) + amino-acid property change
    - Severity classification (conservative / moderate / radical)

This module exposes a small set of clean, stable functions that the
Streamlit front-end (app.py, owned by Member 2 - Mainyatha Paladugu)
imports directly. Nothing in here touches the UI.
"""

from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional

from Bio import Align
from Bio.Align import substitution_matrices


# ---------------------------------------------------------------------------
# Amino-acid property groups.
# Used to flag whether a substitution crosses into a different physicochemical
# class (e.g. acidic -> nonpolar), independent of the raw BLOSUM62 score.
# ---------------------------------------------------------------------------
AA_PROPERTIES: Dict[str, str] = {
    # Nonpolar / hydrophobic
    "G": "nonpolar", "A": "nonpolar", "V": "nonpolar", "L": "nonpolar",
    "I": "nonpolar", "P": "nonpolar", "F": "nonpolar", "M": "nonpolar",
    "W": "nonpolar",
    # Polar, uncharged
    "S": "polar", "T": "polar", "C": "polar", "Y": "polar",
    "N": "polar", "Q": "polar",
    # Positively charged (basic)
    "K": "basic", "R": "basic", "H": "basic",
    # Negatively charged (acidic)
    "D": "acidic", "E": "acidic",
}

SEVERITY_LEVELS = ["conservative", "moderate", "radical"]


@dataclass
class Mutation:
    """A single detected difference between the two aligned sequences."""
    position: int                      # 1-based position, numbered along seq1 (reference)
    original: str                      # residue in seq1 ("-" if this is an insertion)
    mutated: str                       # residue in seq2 ("-" if this is a deletion)
    kind: str                          # "substitution" | "insertion" | "deletion"
    blosum_score: Optional[float] = None
    property_change: Optional[str] = None
    severity: Optional[str] = None

    def to_dict(self) -> Dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# 1. Alignment
# ---------------------------------------------------------------------------
def align_sequences(
    seq1: str,
    seq2: str,
    matrix_name: str = "BLOSUM62",
    open_gap_score: float = -10.0,
    extend_gap_score: float = -0.5,
) -> Tuple[str, str, float]:
    """
    Global pairwise alignment (Needleman-Wunsch) of two protein sequences,
    scored with a substitution matrix (default BLOSUM62).

    Returns (aligned_seq1, aligned_seq2, alignment_score), where the aligned
    sequences are the same length and use "-" for gaps.
    """
    seq1 = _clean_sequence(seq1)
    seq2 = _clean_sequence(seq2)
    if not seq1 or not seq2:
        raise ValueError("Both sequences must be non-empty protein sequences.")

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    alignments = aligner.align(seq1, seq2)
    best = alignments[0]
    aligned1, aligned2 = str(best[0]), str(best[1])
    return aligned1, aligned2, float(best.score)


def _clean_sequence(seq: str) -> str:
    """Strip whitespace/FASTA headers and uppercase a raw sequence string."""
    lines = [line.strip() for line in seq.strip().splitlines()]
    lines = [line for line in lines if line and not line.startswith(">")]
    return "".join(lines).upper()


# ---------------------------------------------------------------------------
# 2. Mutation detection
# ---------------------------------------------------------------------------
def detect_mutations(aligned1: str, aligned2: str) -> List[Mutation]:
    """
    Walk an alignment (two equal-length gapped strings) and record every
    position where the two sequences disagree, as a Mutation.
    """
    if len(aligned1) != len(aligned2):
        raise ValueError("Aligned sequences must be the same length.")

    mutations: List[Mutation] = []
    ref_pos = 0  # 1-based position along seq1 (the reference / wild type)

    for a, b in zip(aligned1, aligned2):
        if a != "-":
            ref_pos += 1
        if a == b:
            continue
        if a == "-":
            mutations.append(Mutation(position=ref_pos, original="-", mutated=b, kind="insertion"))
        elif b == "-":
            mutations.append(Mutation(position=ref_pos, original=a, mutated="-", kind="deletion"))
        else:
            mutations.append(Mutation(position=ref_pos, original=a, mutated=b, kind="substitution"))
    return mutations


# ---------------------------------------------------------------------------
# 3. Impact scoring
# ---------------------------------------------------------------------------
def classify_severity(score: float, property_changed: bool) -> str:
    """
    Turn a raw BLOSUM62 score (+ whether the amino-acid property class
    changed) into a 3-tier severity label.
    """
    if score >= 1:
        level = 0        # conservative
    elif score >= -2:
        level = 1        # moderate
    else:
        level = 2        # radical

    if property_changed:
        level = min(level + 1, len(SEVERITY_LEVELS) - 1)

    return SEVERITY_LEVELS[level]


def score_mutation(mutation: Mutation, matrix_name: str = "BLOSUM62") -> Mutation:
    """Fill in blosum_score / property_change / severity for one mutation."""
    if mutation.kind != "substitution":
        mutation.blosum_score = None
        mutation.property_change = "n/a (indel)"
        mutation.severity = "structural"
        return mutation

    matrix = substitution_matrices.load(matrix_name)
    try:
        score = matrix[mutation.original, mutation.mutated]
    except KeyError:
        score = matrix[mutation.mutated, mutation.original]
    mutation.blosum_score = float(score)

    prop_a = AA_PROPERTIES.get(mutation.original, "unknown")
    prop_b = AA_PROPERTIES.get(mutation.mutated, "unknown")
    mutation.property_change = "unchanged" if prop_a == prop_b else f"{prop_a} -> {prop_b}"

    mutation.severity = classify_severity(mutation.blosum_score, prop_a != prop_b)
    return mutation


def score_mutations(
    mutations: List[Mutation], matrix_name: str = "BLOSUM62"
) -> Tuple[List[Mutation], float, float]:
    """
    Score every mutation in place and combine substitution scores into an
    overall mutation-impact score for the sequence pair.

    Returns (scored_mutations, total_impact_score, average_impact_score).
    Indels are reported but excluded from the numeric impact score, since
    BLOSUM62 has no meaningful entry for a gap.
    """
    scored = [score_mutation(m, matrix_name) for m in mutations]
    substitutions = [m for m in scored if m.kind == "substitution"]

    total = sum(m.blosum_score for m in substitutions) if substitutions else 0.0
    average = total / len(substitutions) if substitutions else 0.0
    return scored, total, average


# ---------------------------------------------------------------------------
# 4. Convenience: one-call entry point for the UI
# ---------------------------------------------------------------------------
def analyze(seq1: str, seq2: str, matrix_name: str = "BLOSUM62") -> Dict:
    """
    Run the full pipeline: align -> detect mutations -> score -> summarize.

    This is the single function Partner B's Streamlit app needs to call.
    """
    aligned1, aligned2, alignment_score = align_sequences(seq1, seq2, matrix_name)
    mutations = detect_mutations(aligned1, aligned2)
    scored_mutations, total_impact, average_impact = score_mutations(mutations, matrix_name)

    n_sub = sum(1 for m in scored_mutations if m.kind == "substitution")
    n_ins = sum(1 for m in scored_mutations if m.kind == "insertion")
    n_del = sum(1 for m in scored_mutations if m.kind == "deletion")

    return {
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "alignment_score": alignment_score,
        "mutations": scored_mutations,
        "mutations_as_dicts": [m.to_dict() for m in scored_mutations],
        "num_mutations": len(scored_mutations),
        "num_substitutions": n_sub,
        "num_insertions": n_ins,
        "num_deletions": n_del,
        "total_impact_score": total_impact,
        "average_impact_score": average_impact,
    }


if __name__ == "__main__":
    # Quick manual smoke test (see tests/test_mutation_scorer.py for real tests)
    seq1 = "MVLSPADKTNVKAAWGKVGAH"
    seq2 = "MVLSPADKTNVKAAWGEVGAH"  # one substitution: A -> E
    result = analyze(seq1, seq2)
    print("Seq1:", result["aligned_seq1"])
    print("Seq2:", result["aligned_seq2"])
    print(f"Alignment score: {result['alignment_score']}")
    for m in result["mutations"]:
        print(f"  pos {m.position}: {m.original}->{m.mutated} "
              f"({m.kind}, BLOSUM62={m.blosum_score}, {m.property_change}, {m.severity})")
    print(f"Total impact score: {result['total_impact_score']}")
    print(f"Average impact score: {result['average_impact_score']}")
