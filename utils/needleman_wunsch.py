from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal, Tuple


SequenceType = Literal["DNA", "RNA", "Protein"]


DNA_ALPHABET = set("ATCG")
RNA_ALPHABET = set("AUCG")
# 20 standard amino acids (no B, J, O, U, X, Z)
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

# IUPAC ambiguity codes for nucleotides (DNA/RNA)
NUCLEOTIDE_AMBIGUOUS = set("RYSWKMBDHVN")
# Common protein ambiguity codes (B = D/N, Z = E/Q, X = unknown)
PROTEIN_AMBIGUOUS = set("BXZ")


@dataclass(frozen=True)
class AlignmentResult:
    aligned_seq1: str
    aligned_seq2: str
    score: int
    matrix: List[List[int]]
    traceback_path: List[Tuple[int, int]]  # (i, j) coordinates in DP matrix
    match_line: str
    matches: int
    mismatches: int
    gaps: int
    identity_pct: float


def _strip_fasta(seq: str) -> str:
    """
    Accepts raw user input and removes FASTA headers (lines starting with '>').
    """
    lines = str(seq).splitlines()
    kept = [ln for ln in lines if not ln.lstrip().startswith(">")]
    return "\n".join(kept)


def _clean_sequence(seq: str) -> str:
    # Remove FASTA headers, then remove all whitespace/newlines, then uppercase
    no_headers = _strip_fasta(seq)
    return "".join(no_headers.split()).upper()


def validate_sequence(
    seq: str,
    sequence_type: SequenceType,
    *,
    allow_ambiguous: bool = False,
) -> Tuple[bool, str | None, str]:
    """
    Returns: (is_valid, error_message, cleaned_seq)
    """
    cleaned = _clean_sequence(seq)
    if not cleaned:
        return False, "Sequence cannot be empty.", cleaned

    if sequence_type == "DNA":
        alphabet = set(DNA_ALPHABET) | (NUCLEOTIDE_AMBIGUOUS if allow_ambiguous else set())
    elif sequence_type == "RNA":
        alphabet = set(RNA_ALPHABET) | (NUCLEOTIDE_AMBIGUOUS if allow_ambiguous else set())
    else:
        alphabet = set(PROTEIN_ALPHABET) | (PROTEIN_AMBIGUOUS if allow_ambiguous else set())
    invalid = sorted({c for c in cleaned if c not in alphabet})
    if invalid:
        allowed = "".join(sorted(alphabet))
        return (
            False,
            f"Invalid character(s): {', '.join(invalid)}. Allowed for {sequence_type}: {allowed}.",
            cleaned,
        )
    return True, None, cleaned


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> AlignmentResult:
    """
    Global alignment (Needleman–Wunsch) with constant gap penalty.

    Scoring:
      - match: +match_score
      - mismatch: +mismatch_penalty (typically negative)
      - gap: +gap_penalty (typically negative)
    """
    s1 = _clean_sequence(seq1)
    s2 = _clean_sequence(seq2)

    n = len(s1)
    m = len(s2)

    # DP matrix size (n+1) x (m+1)
    F: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]
    # Trace matrix: store direction(s) chosen; use ints for compactness
    # 0 = diag, 1 = up, 2 = left
    T: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    # init
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + gap_penalty
        T[i][0] = 1  # up
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + gap_penalty
        T[0][j] = 2  # left

    # fill
    for i in range(1, n + 1):
        c1 = s1[i - 1]
        for j in range(1, m + 1):
            c2 = s2[j - 1]
            diag = F[i - 1][j - 1] + (match_score if c1 == c2 else mismatch_penalty)
            up = F[i - 1][j] + gap_penalty
            left = F[i][j - 1] + gap_penalty

            best = diag
            direction = 0
            if up > best:
                best = up
                direction = 1
            if left > best:
                best = left
                direction = 2
            # tie-breaker: prefer diag, then up, then left (keeps alignments intuitive)
            if best == diag:
                direction = 0
            elif best == up:
                direction = 1
            else:
                direction = 2

            F[i][j] = best
            T[i][j] = direction

    # traceback
    i, j = n, m
    a1: List[str] = []
    a2: List[str] = []
    path: List[Tuple[int, int]] = [(i, j)]

    while i > 0 or j > 0:
        if i > 0 and j > 0 and T[i][j] == 0:
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or T[i][j] == 1):
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        else:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1
        path.append((i, j))

    a1.reverse()
    a2.reverse()
    path.reverse()  # from (0,0) to (n,m)

    aligned1 = "".join(a1)
    aligned2 = "".join(a2)

    # metrics + match-line
    match_line_chars: List[str] = []
    matches = 0
    mismatches = 0
    gaps = 0

    for c1, c2 in zip(aligned1, aligned2):
        if c1 == "-" or c2 == "-":
            gaps += 1
            match_line_chars.append(" ")
        elif c1 == c2:
            matches += 1
            match_line_chars.append("|")
        else:
            mismatches += 1
            match_line_chars.append(" ")

    aligned_len = max(1, len(aligned1))
    identity_pct = (matches / aligned_len) * 100.0

    return AlignmentResult(
        aligned_seq1=aligned1,
        aligned_seq2=aligned2,
        score=F[n][m],
        matrix=F,
        traceback_path=path,
        match_line="".join(match_line_chars),
        matches=matches,
        mismatches=mismatches,
        gaps=gaps,
        identity_pct=identity_pct,
    )

