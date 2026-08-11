#!/usr/bin/env python3
"""
Truncate BAM reads to contain ONLY the bases that align within the Gviz
plotting windows defined below.

Unlike simple soft-clipping, this physically shortens the query sequence so
the read has no bases outside the window at all.  This is required because
Gviz's sequenceLayer iterates over the full query sequence (including
soft-clipped bases) when computing mismatch positions, causing a
"replaceAt: ranges off-limits" crash for reads longer than the tiny (35-42
bp) IDH1/IDH2/TERTp plotting windows.

Clinical samples are unaffected (they have zero coverage at these windows).
This script is specific to PAU59949 (WGS cell-line, 40-70 kb reads).

Usage:
    python3 clip_bam_to_plot_windows.py <input.bam> <output.bam>
"""
import pysam
import sys


# BAM reference_start is 0-based; R plotTracks(from, to) is 1-based inclusive.
# Conversion: 0-based_start = R_from - 1,  0-based_exclusive_end = R_to
#   (the exclusive end in 0-based equals the inclusive end in 1-based numerically)
WINDOWS = [
    ("chr2",  208248369, 208248405),   # IDH1: R from=208248370, to=208248405
    ("chr15", 90088586,  90088622),    # IDH2: R from=90088587,  to=90088622
    ("chr5",  1295102,   1295145),     # TERT: R from=1295103,   to=1295145
    ("chr7",  55019016,  55211628),    # EGFR: R from=55019017,  to=55211628
]


def find_window(chrom, ref_start, ref_end):
    for wc, ws, we in WINDOWS:
        if chrom == wc and ref_start < we and ref_end > ws:
            return ws, we
    return None


def truncate_read_to_window(read, window_start, window_end):
    """
    Return a modified read whose query sequence contains only the bases that
    align to reference positions within [window_start, window_end).

    The resulting read has:
      - A pure M cigar (no soft clips, no indel ops)
      - query_sequence and query_qualities trimmed to those bases only
      - reference_start set to the first in-window reference position
      - MD tag removed (would be inconsistent with the new cigar)

    Insertions within the window consume extra query bases but no reference
    bases; we keep only the first `middle_r` query bases to avoid exceeding
    the reference span.  Deletions (skipped reference positions) are dropped
    from the cigar for simplicity — acceptable for a display-only BAM.

    Returns None if the read has no bases overlapping the window.
    """
    if read.is_unmapped or not read.cigartuples or not read.query_sequence:
        return read

    seq  = read.query_sequence
    qual = read.query_qualities          # numpy array or None

    ref_start = read.reference_start     # 0-based
    ref_end   = read.reference_end       # 0-based exclusive (may be None)

    # ---- find first / last query positions that map to in-window ref pos ----
    first_q     = None
    last_q      = None
    new_r_start = None

    for qpos, rpos in read.get_aligned_pairs():
        if rpos is None:
            # Insertion: no reference position.
            # Keep if we are already inside the window.
            if first_q is not None and qpos is not None:
                last_q = qpos
            continue
        if qpos is None:
            # Deletion: no query base — skip.
            continue
        if window_start <= rpos < window_end:
            if first_q is None:
                first_q     = qpos
                new_r_start = rpos
            last_q = qpos

    if first_q is None:
        return None   # no aligned bases in this window

    # ---- how many reference positions does this window span? ---------------
    # Use new_r_start (the actual first in-window position of THIS read), not
    # window_start.  A read that begins inside the window (e.g. at window+1)
    # must be limited to window_end - new_r_start bases, not window_end -
    # window_start, otherwise the M cigar overshoots the window by 1 and Gviz
    # computes an out-of-bounds relative position (window_end + 1 - from + 1).
    clipped_ref_end = min(ref_end, window_end) if ref_end else window_end
    middle_r        = clipped_ref_end - new_r_start   # ref bases from actual start to window end
    middle_q        = last_q - first_q + 1            # query bases selected

    # Keep at most middle_r query bases so the M cigar does not overshoot
    # the window on the right (handles insertions by trimming them away).
    used_q = min(middle_q, middle_r)
    if used_q <= 0:
        return None

    new_seq  = seq[first_q : first_q + used_q]
    new_qual = qual[first_q : first_q + used_q] if qual is not None else None

    # ---- rewrite the read fields -----------------------------------------
    # Order matters in pysam: set cigar/pos first, then seq, then qual.
    read.cigartuples      = [(0, used_q)]   # pure M — no soft clips
    read.reference_start  = new_r_start
    read.query_sequence   = new_seq         # physically shorter sequence
    if new_qual is not None:
        read.query_qualities = new_qual

    # Remove MD tag — it described the original full read alignment.
    if read.has_tag("MD"):
        read.set_tag("MD", None)

    return read


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bam> <output.bam>", file=sys.stderr)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]
    tmp_path = out_path + ".unsorted.tmp.bam"
    n_truncated = 0
    n_unchanged = 0

    with pysam.AlignmentFile(in_path, "rb") as bam_in, \
         pysam.AlignmentFile(tmp_path, "wb", header=bam_in.header) as bam_out:

        for read in bam_in:
            if read.is_unmapped or not read.cigartuples or not read.reference_name:
                bam_out.write(read)
                n_unchanged += 1
                continue

            window = find_window(
                read.reference_name,
                read.reference_start,
                read.reference_end,
            )

            if window is None:
                bam_out.write(read)
                n_unchanged += 1
                continue

            truncated = truncate_read_to_window(read, window[0], window[1])
            if truncated is not None:
                bam_out.write(truncated)
                n_truncated += 1
            # reads with zero in-window bases are dropped

    print(f"Truncated: {n_truncated}  Unchanged: {n_unchanged}", file=sys.stderr)

    # Sort and index using pysam's bundled samtools
    pysam.sort("-o", out_path, tmp_path)
    pysam.index(out_path)

    import os
    os.unlink(tmp_path)
    print(f"Sorted and indexed: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
