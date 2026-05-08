from __future__ import annotations

from flask import Flask, render_template, request

from utils.needleman_wunsch import (
    needleman_wunsch,
    validate_sequence,
)

app = Flask(__name__)


def _parse_int(value: str | None, default: int) -> int:
    if value is None or str(value).strip() == "":
        return default
    return int(value)


@app.get("/")
def index():
    return render_template("index.html")


@app.get("/about")
def about():
    return render_template("about.html")


@app.post("/analyze")
def analyze():
    seq1_raw = request.form.get("sequence1", "")
    seq2_raw = request.form.get("sequence2", "")
    seq_type = request.form.get("sequence_type", "DNA")
    allow_ambiguous = request.form.get("allow_ambiguous") == "on"

    errors: list[str] = []

    try:
        match_score = _parse_int(request.form.get("match_score"), 1)
        mismatch_penalty = _parse_int(request.form.get("mismatch_penalty"), -1)
        gap_penalty = _parse_int(request.form.get("gap_penalty"), -2)
    except ValueError:
        errors.append("Scoring parameters must be integers (e.g., match=1, mismatch=-1, gap=-2).")
        match_score, mismatch_penalty, gap_penalty = 1, -1, -2

    if seq_type not in ("DNA", "RNA", "Protein"):
        errors.append("Sequence type is invalid.")
        seq_type = "DNA"

    ok1, err1, seq1 = validate_sequence(seq1_raw, seq_type, allow_ambiguous=allow_ambiguous)  # type: ignore[arg-type]
    ok2, err2, seq2 = validate_sequence(seq2_raw, seq_type, allow_ambiguous=allow_ambiguous)  # type: ignore[arg-type]
    if not ok1 and err1:
        errors.append(f"Sequence 1: {err1}")
    if not ok2 and err2:
        errors.append(f"Sequence 2: {err2}")

    if errors:
        return render_template(
            "index.html",
            errors=errors,
            form_data={
                "sequence1": seq1_raw,
                "sequence2": seq2_raw,
                "sequence_type": seq_type,
                "allow_ambiguous": allow_ambiguous,
                "match_score": request.form.get("match_score", "1"),
                "mismatch_penalty": request.form.get("mismatch_penalty", "-1"),
                "gap_penalty": request.form.get("gap_penalty", "-2"),
            },
        )

    result = needleman_wunsch(
        seq1=seq1,
        seq2=seq2,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
    )

    # Build axis labels for matrix display
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)
    path_set = set(result.traceback_path)

    return render_template(
        "results.html",
        seq_type=seq_type,
        input_seq1=seq1,
        input_seq2=seq2,
        scoring={
            "match_score": match_score,
            "mismatch_penalty": mismatch_penalty,
            "gap_penalty": gap_penalty,
        },
        alignment={
            "aligned1": result.aligned_seq1,
            "match_line": result.match_line,
            "aligned2": result.aligned_seq2,
            "score": result.score,
        },
        metrics={
            "matches": result.matches,
            "mismatches": result.mismatches,
            "gaps": result.gaps,
            "identity_pct": result.identity_pct,
            "aligned_length": len(result.aligned_seq1),
        },
        matrix=result.matrix,
        row_labels=row_labels,
        col_labels=col_labels,
        traceback_path=result.traceback_path,
        traceback_set=path_set,
        final_cell=(len(seq1), len(seq2)),
    )


if __name__ == "__main__":
    app.run(debug=True)

