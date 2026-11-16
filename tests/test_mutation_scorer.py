"""
tests/test_mutation_scorer.py

Unit tests for the core bioinformatics engine, using a mix of synthetic
sequences (to pin down exact behaviour) and a real, published protein
variant (to sanity-check the pipeline biologically).

Run with:  pytest tests/ -v
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mutation_scorer import (
    align_sequences,
    detect_mutations,
    score_mutations,
    classify_severity,
    analyze,
)


# ---------------------------------------------------------------------------
# Alignment + mutation detection
# ---------------------------------------------------------------------------
def test_identical_sequences_have_no_mutations():
    seq = "MVLSPADKTNVKAAWGKVGAH"
    result = analyze(seq, seq)
    assert result["num_mutations"] == 0
    assert result["total_impact_score"] == 0.0


def test_single_point_substitution_is_detected():
    # Example from the team's own scratch script: one substitution (K17E)
    seq1 = "MVLSPADKTNVKAAWGKVGAH"
    seq2 = "MVLSPADKTNVKAAWGEVGAH"
    result = analyze(seq1, seq2)

    assert result["num_mutations"] == 1
    assert result["num_substitutions"] == 1

    mutation = result["mutations"][0]
    assert mutation.kind == "substitution"
    assert mutation.original == "K"
    assert mutation.mutated == "E"
    assert mutation.position == 17
    assert mutation.blosum_score is not None


def test_insertion_is_detected():
    seq1 = "ACDEFG"
    seq2 = "ACDWEFG"  # W inserted between D and E
    aligned1, aligned2, _ = align_sequences(seq1, seq2)
    mutations = detect_mutations(aligned1, aligned2)

    assert len(mutations) == 1
    assert mutations[0].kind == "insertion"
    assert mutations[0].mutated == "W"


def test_deletion_is_detected():
    seq1 = "ACDWEFG"
    seq2 = "ACDEFG"  # W deleted relative to seq1
    aligned1, aligned2, _ = align_sequences(seq1, seq2)
    mutations = detect_mutations(aligned1, aligned2)

    assert len(mutations) == 1
    assert mutations[0].kind == "deletion"
    assert mutations[0].original == "W"


# ---------------------------------------------------------------------------
# Impact scoring / severity classification
# ---------------------------------------------------------------------------
def test_classify_severity_thresholds():
    assert classify_severity(score=4, property_changed=False) == "conservative"
    assert classify_severity(score=0, property_changed=False) == "moderate"
    assert classify_severity(score=-5, property_changed=False) == "radical"
    # a property-class change escalates severity by one tier
    assert classify_severity(score=4, property_changed=True) == "moderate"


def test_real_world_sickle_cell_mutation():
    """
    Human hemoglobin subunit beta (HBB), UniProt P68871.
    The sickle-cell variant is the classic Glu6Val (E6V) substitution
    (numbered from the mature chain; position 7 here because the
    UniProt sequence still includes the initiator Met at position 1).
    """
    wild_type = (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV"
        "KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKE"
        "FTPPVQAAYQKVVAGVANALAHKYH"
    )
    sickle_variant = wild_type[:6] + "V" + wild_type[7:]  # Glu -> Val at position 7

    result = analyze(wild_type, sickle_variant)

    assert result["num_mutations"] == 1
    mutation = result["mutations"][0]
    assert mutation.kind == "substitution"
    assert mutation.position == 7
    assert mutation.original == "E"
    assert mutation.mutated == "V"
    # acidic -> nonpolar is a real physicochemical class change
    assert mutation.property_change == "acidic -> nonpolar"
    # BLOSUM62 scores E->V negatively, and the property class changes,
    # so this should come out as a high-severity ("radical") change --
    # consistent with sickle-cell disease being caused by exactly this
    # substitution replacing a surface-exposed, charged residue with a
    # hydrophobic one.
    assert mutation.severity == "radical"


def test_total_impact_score_sums_substitutions_only():
    seq1 = "ACDWEFG"
    seq2 = "ACDEFH"  # 1 deletion (W) + 1 substitution (G->H); lengths differ
    aligned1, aligned2, _ = align_sequences(seq1, seq2)
    mutations = detect_mutations(aligned1, aligned2)
    scored, total, average = score_mutations(mutations)

    substitutions = [m for m in scored if m.kind == "substitution"]
    indels = [m for m in scored if m.kind != "substitution"]

    assert len(substitutions) >= 1
    assert all(m.blosum_score is None for m in indels)
    assert total == sum(m.blosum_score for m in substitutions)
    if substitutions:
        assert average == total / len(substitutions)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))
