"""
app.py
------
Streamlit website for the Protein Sequence Mutation Scorer.

What it does:
  1. Lets the user paste a normal (reference) protein sequence and a changed one.
  2. Calls analyze() from core/mutation_scorer.py to align them and score the mutations.
  3. Shows the aligned sequences (differences highlighted), a table of mutations,
     and a bar chart of the BLOSUM62 scores.
"""

import pandas as pd
import streamlit as st

from core.mutation_scorer import analyze

VALID_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")

# Example: human hemoglobin beta (HBB) normal vs. sickle-cell variant (Glu -> Val)
EXAMPLE_NORMAL = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV"
    "KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKE"
    "FTPPVQAAYQKVVAGVANALAHKYH"
)
EXAMPLE_CHANGED = EXAMPLE_NORMAL[:6] + "V" + EXAMPLE_NORMAL[7:]

SEVERITY_COLORS = {
    "conservative": "#2e7d32",
    "moderate": "#ef6c00",
    "radical": "#c62828",
    "structural": "#6a1b9a",
}


# ---------------------------------------------------------------------------
# Small helper functions
# ---------------------------------------------------------------------------
def load_example():
    """Fill both text boxes with the sickle-cell example."""
    st.session_state["seq1"] = EXAMPLE_NORMAL
    st.session_state["seq2"] = EXAMPLE_CHANGED


def clean(text):
    """Remove FASTA header lines (starting with >), spaces and line breaks."""
    lines = [ln.strip() for ln in text.splitlines() if not ln.startswith(">")]
    return "".join(lines).upper().replace(" ", "")


def highlight_alignment(aligned1, aligned2, width=60):
    """Return HTML showing the two aligned sequences, with differences in red."""
    blocks = []
    for start in range(0, len(aligned1), width):
        a = aligned1[start:start + width]
        b = aligned2[start:start + width]
        top, bottom, marks = "", "", ""
        for x, y in zip(a, b):
            if x == y:
                top += x
                bottom += y
                marks += "&nbsp;"
            else:
                top += f"<span style='background:#ffcdd2;color:#b71c1c'>{x}</span>"
                bottom += f"<span style='background:#ffcdd2;color:#b71c1c'>{y}</span>"
                marks += "*"
        blocks.append(
            f"<div style='margin-bottom:12px'>"
            f"<div>Normal &nbsp;{top}</div>"
            f"<div>Changed {bottom}</div>"
            f"<div style='color:#b71c1c'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{marks}</div>"
            f"</div>"
        )
    return "<div style='font-family:monospace;font-size:14px;overflow-x:auto'>" + "".join(blocks) + "</div>"


# ---------------------------------------------------------------------------
# Page layout
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Protein Mutation Scorer", page_icon="🧬", layout="wide")

st.title("🧬 Protein Sequence Mutation Scorer")
st.write(
    "Paste a normal protein sequence and a changed (mutated) one. The tool lines them up, "
    "finds every difference, and estimates how serious each change might be."
)

st.button("Load example (sickle-cell hemoglobin)", on_click=load_example)

col1, col2 = st.columns(2)
with col1:
    seq1 = st.text_area("Normal (reference) sequence", key="seq1", height=180,
                        placeholder="Paste amino acid letters here, e.g. MVLSPADKTN...")
with col2:
    seq2 = st.text_area("Changed (mutated) sequence", key="seq2", height=180,
                        placeholder="Paste the changed sequence here")

if st.button("Analyze", type="primary"):
    s1, s2 = clean(seq1), clean(seq2)

    # --- Check the input before running anything ---
    if not s1 or not s2:
        st.warning("Please enter both sequences.")
        st.stop()
    bad = sorted((set(s1) | set(s2)) - VALID_LETTERS)
    if bad:
        st.error(f"These characters are not standard amino acid letters: {', '.join(bad)}")
        st.stop()

    # --- Run the main analysis from core/mutation_scorer.py ---
    result = analyze(s1, s2)

    # --- Summary numbers ---
    st.subheader("Summary")
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Total differences", result["num_mutations"])
    c2.metric("Substitutions", result["num_substitutions"])
    c3.metric("Insertions", result["num_insertions"])
    c4.metric("Deletions", result["num_deletions"])
    c5.metric("Average BLOSUM62 score", f"{result['average_impact_score']:.2f}")
    st.caption(
        "BLOSUM62 score: higher = the swap is common in nature (probably harmless); "
        "lower or negative = the swap is rare (more likely to matter)."
    )

    if result["num_mutations"] == 0:
        st.success("The two sequences are identical. No mutations found.")
        st.stop()

    # --- Aligned sequences with differences highlighted ---
    st.subheader("Aligned sequences")
    st.markdown(
        highlight_alignment(result["aligned_seq1"], result["aligned_seq2"]),
        unsafe_allow_html=True,
    )
    st.caption("Red letters and * mark the differences. A '-' is a gap (an insertion or deletion).")

    # --- Table of mutations ---
    st.subheader("Mutations found")
    table = pd.DataFrame(result["mutations_as_dicts"])
    table = table.rename(columns={
        "position": "Position",
        "original": "Original",
        "mutated": "Changed to",
        "kind": "Type",
        "blosum_score": "BLOSUM62 score",
        "property_change": "Property change",
        "severity": "Severity",
    })
    st.dataframe(table, width="stretch", hide_index=True)

    # --- Chart of scores (substitutions only) ---
    subs = [m for m in result["mutations"] if m.kind == "substitution"]
    if subs:
        st.subheader("BLOSUM62 score of each substitution")
        chart_data = pd.DataFrame(
            {"BLOSUM62 score": [m.blosum_score for m in subs]},
            index=[f"{m.original}{m.position}{m.mutated}" for m in subs],
        )
        st.bar_chart(chart_data)

    # --- Severity labels in plain words ---
    st.subheader("Severity")
    for m in result["mutations"]:
        color = SEVERITY_COLORS.get(m.severity, "#444")
        name = f"{m.original}{m.position}{m.mutated}" if m.kind == "substitution" else (
            f"{m.kind} at position {m.position} ({m.original}{m.mutated})".replace("-", "")
        )
        st.markdown(
            f"**{name}** — <span style='color:{color};font-weight:bold'>{m.severity}</span> "
            f"({m.property_change})",
            unsafe_allow_html=True,
        )
