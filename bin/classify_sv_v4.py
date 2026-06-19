#!/usr/bin/env python3
"""
classify_sv_v4.py — Classify Sniffles SVs by gene-level impact using BEDTools.

v2: all-genes intrachromosomal fusion detection (CHIC2 inside KIT fix).
v3: --bp-window to expand breakpoint intervals in fusion detection.
v4: adds gene_type_A / gene_type_B columns to sv_fusions.tsv and a new
    output sv_fusions_gbm_protein.tsv — fusions where at least one gene is
    in the GBM driver list AND both partners are protein_coding.

Priority hierarchy:
  full_gene_loss > exon_disruption > intronic_disruption > intergenic

BND records (translocations) are skipped from interval intersections and
reported separately as 'unclassified_bnd'.

Usage:
    python classify_sv.py \\
        --vcf sample.sniffles.vcf.gz \\
        --gff3 gencode.v48.annotation.gff3 \\
        --out sv_classified.tsv \\
        [--gene-list gbm_genes.txt --out-genes sv_gbm.tsv] \\
        [--protein-coding-only]

--gene-list   Text file with gene names (one per line, optional line-number prefix).
              When provided, a second TSV is written (--out-genes) containing only SVs
              that overlap at least one gene in the list. An extra column 'gbm_genes'
              shows which listed genes are hit.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def run(cmd, check=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and r.returncode != 0:
        sys.exit(f"Command failed:\n  {cmd}\nstderr:\n  {r.stderr.strip()}")
    return r.stdout


def extract_sv_bed(vcf, sv_bed, bnd_bed):
    """Convert VCF to BED (0-based). BNDs written to separate file."""
    raw = run(
        f"bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%ID\\t%INFO/SVTYPE\\n' {vcf}"
    )
    with open(sv_bed, "w") as sf, open(bnd_bed, "w") as bf:
        for line in raw.splitlines():
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            chrom, pos, end, sv_id, svtype = cols
            if svtype == "BND" or end == ".":
                bf.write(f"{chrom}\t{int(pos) - 1}\t{pos}\t{sv_id}\t{svtype}\n")
                continue
            try:
                start0 = int(pos) - 1
                end_i = int(end)
                if end_i <= start0:  # point event (INS with END=POS)
                    end_i = start0 + 1
                sf.write(f"{chrom}\t{start0}\t{end_i}\t{sv_id}\t{svtype}\n")
            except ValueError:
                continue


def _vcf_has_info_tag(vcf, tag):
    """Check the VCF header for an INFO/<tag> definition."""
    r = subprocess.run(f"bcftools view -h {vcf}", shell=True,
                       capture_output=True, text=True)
    return bool(re.search(rf"^##INFO=<ID={re.escape(tag)},", r.stdout, re.MULTILINE))


def extract_sv_info(vcf):
    """Return {sv_id: (support, af)} keyed by the Sniffles2 record ID
    (e.g. Sniffles2.INV.xxx), read from INFO/SUPPORT and INFO/AF.

    Querying a tag absent from the VCF header makes bcftools query fail
    outright (empty stdout), which would silently drop SUPPORT too — so
    only query AF when it's actually declared in the header, otherwise
    leave it as ".".
    """
    af_tag = "%INFO/AF" if _vcf_has_info_tag(vcf, "AF") else "."

    raw = run(
        f"bcftools query -f '%ID\\t%INFO/SUPPORT\\t{af_tag}\\n' {vcf}",
        check=False,
    )
    info = {}
    for line in raw.splitlines():
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        sv_id, support, af = cols[0], cols[1], cols[2]
        info[sv_id] = (support, af)
    return info


def extract_annotation_beds(gff3, genes_bed, exons_bed, exons_detail_bed,
                            gene_types=None):
    """Parse GFF3 → genes.bed, exons.bed, and exons_detail.bed.

    exons_detail_bed name col = 'gene_name|transcript_id|exon_number'
    (used for per-breakpoint annotation).
    """
    with open(gff3) as src, \
         open(genes_bed, "w") as gf, \
         open(exons_bed, "w") as ef, \
         open(exons_detail_bed, "w") as ed:
        for line in src:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = cols
            if feat not in ("gene", "exon"):
                continue
            attr = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr[k] = v
            gtype = attr.get("gene_type", "unknown")
            if gene_types and gtype not in gene_types:
                continue
            gene_name = attr.get("gene_name", attr.get("gene_id", "unknown"))
            gene_id = attr.get("gene_id", "unknown")
            name = f"{gene_name}|{gene_id}|{gtype}"
            start0 = int(start) - 1
            bed_line = f"{chrom}\t{start0}\t{end}\t{name}\t.\t{strand}\n"
            if feat == "gene":
                gf.write(bed_line)
            else:
                ef.write(bed_line)
                # detailed exon BED: name = gene_name|transcript_id|exon_number
                transcript_id = attr.get("transcript_id", "unknown")
                exon_num = attr.get("exon_number", "?")
                # include gene_type so annotate_positions can rank by biotype
                detail_name = f"{gene_name}|{transcript_id}|{exon_num}|{gtype}"
                ed.write(f"{chrom}\t{start0}\t{end}\t{detail_name}\t.\t{strand}\n")


def intersect(a, b, extra=""):
    """Run bedtools intersect -wa -wb; return list of split-column lines."""
    cmd = f"bedtools intersect -a {a} -b {b} {extra} -wa -wb"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode not in (0, 1):  # 1 = no overlaps found
        sys.exit(f"bedtools intersect failed:\n  {r.stderr.strip()}")
    return [ln.split("\t") for ln in r.stdout.splitlines() if ln]


def gene_name_from_field(field):
    return field.split("|")[0]


def build_map_sv_a(lines, sv_id_col, gene_col):
    """Build sv_id → set(gene_name) when SV is the -a file."""
    m = defaultdict(set)
    for cols in lines:
        if len(cols) > max(sv_id_col, gene_col):
            m[cols[sv_id_col]].add(gene_name_from_field(cols[gene_col]))
    return m


def build_map_gene_a(lines, sv_id_col, gene_col):
    """Build sv_id → set(gene_name) when gene is the -a file."""
    m = defaultdict(set)
    for cols in lines:
        if len(cols) > max(sv_id_col, gene_col):
            m[cols[sv_id_col]].add(gene_name_from_field(cols[gene_col]))
    return m


def read_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 5:
                rows.append(cols)
    return rows


def make_breakpoints_bed(sv_bed, bp_bed):
    """Write a 1-bp BED for each SV breakpoint (left = POS, right = END).

    Point events (INS, where end = start+1) get only the left breakpoint.
    Name field: sv_id|left or sv_id|right.
    """
    with open(bp_bed, "w") as out:
        for cols in read_bed(sv_bed):
            chrom, start, end, sv_id, svtype = cols[:5]
            s, e = int(start), int(end)
            out.write(f"{chrom}\t{s}\t{s + 1}\t{sv_id}|left\t{svtype}\n")
            if e > s + 1:  # skip right BP for point events
                out.write(f"{chrom}\t{e - 1}\t{e}\t{sv_id}|right\t{svtype}\n")


# Biotype priority for gene selection when multiple genes overlap one position.
# Lower number = higher priority.  protein_coding is always preferred.
_BIOTYPE_PRIORITY = {
    "protein_coding": 0,
    "lncRNA": 2,
    "processed_pseudogene": 3,
    "unprocessed_pseudogene": 3,
}


def _gene_biotype_rank(name_field):
    """Return sort key (priority, gene_name) for a genes.bed name field."""
    parts = name_field.split("|")
    gtype = parts[2] if len(parts) > 2 else "unknown"
    return (_BIOTYPE_PRIORITY.get(gtype, 1), parts[0])


def annotate_positions(pos_bed, exons_detail_bed, genes_bed):
    """Intersect every record in pos_bed with exons and genes.

    Returns {rec_id: (gene_name, element)} where element is
    'exon_N', 'intronic', or 'intergenic'.

    When multiple genes overlap the same position, protein_coding genes are
    preferred over lncRNAs, pseudogenes, and other biotypes so that e.g.
    EP300 is reported rather than EP300-AS1.
    The lowest exon_number across overlapping transcripts is chosen for exons.
    """
    exon_map = defaultdict(list)
    for cols in intersect(pos_bed, exons_detail_bed):
        if len(cols) < 9:
            continue
        rec_id = cols[3]
        parts = cols[8].split("|")
        # parts: gene_name | transcript_id | exon_number | gene_type
        gene_name = parts[0]
        exon_num  = parts[2] if len(parts) > 2 else "?"
        gtype     = parts[3] if len(parts) > 3 else "unknown"
        exon_map[rec_id].append((gene_name, exon_num, gtype))

    # Store full name field so we can rank by biotype
    gene_map = defaultdict(list)
    for cols in intersect(pos_bed, genes_bed):
        if len(cols) < 9:
            continue
        gene_map[cols[3]].append(cols[8])

    result = {}
    for cols in read_bed(pos_bed):
        rec_id = cols[3]
        if rec_id in exon_map:
            # Sort by: biotype priority first, then lowest exon number
            try:
                hits = sorted(
                    exon_map[rec_id],
                    key=lambda x: (_BIOTYPE_PRIORITY.get(x[2], 1), int(x[1]))
                )
            except ValueError:
                hits = sorted(exon_map[rec_id],
                              key=lambda x: _BIOTYPE_PRIORITY.get(x[2], 1))
            gene, exon_num, _ = hits[0]
            result[rec_id] = (gene, f"exon_{exon_num}")
        elif rec_id in gene_map:
            best = min(gene_map[rec_id], key=_gene_biotype_rank)
            result[rec_id] = (gene_name_from_field(best), "intronic")
        else:
            result[rec_id] = ("", "intergenic")
    return result


def annotate_breakpoints(bp_bed, exons_detail_bed, genes_bed, outfile):
    """Write per-breakpoint TSV using annotate_positions."""
    annot = annotate_positions(bp_bed, exons_detail_bed, genes_bed)
    with open(bp_bed) as bf, open(outfile, "w") as out:
        out.write("sv_id\tbp_side\tchrom\tbp_pos\tsvtype\tgene\telement\n")
        for line in bf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, start, _, bp_id, svtype = cols[:5]
            sv_id, side = bp_id.rsplit("|", 1)
            gene, element = annot.get(bp_id, ("", "intergenic"))
            out.write(
                f"{sv_id}\t{side}\t{chrom}\t{start}\t{svtype}\t{gene}\t{element}\n"
            )


def parse_bnd_alt(alt):
    """Return (mate_chrom, mate_pos) from a BND ALT field, or (None, None)."""
    m = re.search(r'[\[\]]([^:]+):(\d+)[\[\]]', alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def detect_fusions(bp_bed, bnd_context, exons_detail_bed, genes_bed,
                   gene_set, outfile, min_fusion_size=10_000, bp_window=0,
                   sv_info=None):
    """Detect gene fusion candidates; write TSV.

    Two sources:
      intrachromosomal — DEL/DUP/INV where the left breakpoint is inside
                         gene A and the right breakpoint is inside gene B != A.
      translocation    — paired BND records where both breakpoints fall
                         inside distinct genes (mates resolved via ALT field).

    bp_window: expand each breakpoint interval by this many bp before
               intersecting with gene annotations (fusion detection only —
               sv_breakpoints.tsv stays exact at 1 bp).

    sv_info: {sv_id: (support, af)} from INFO/SUPPORT and INFO/AF, keyed by
              the Sniffles2 record ID (e.g. Sniffles2.INV.xxx). Looked up by
              sv_id for intrachromosomal fusions and by bnd_id for BNDs.

    Columns: sv_id, fusion_name, svtype, source,
             chrom_A, pos_A, gene_A, element_A, gene_type_A,
             chrom_B, pos_B, gene_B, element_B, gene_type_B,
             both_exonic, both_protein_coding,
             any_in_gene_list, both_in_gene_list,
             support, af
    """
    sv_info = sv_info or {}
    import tempfile as _tf

    # Build gene_name → gene_type lookup from genes_bed
    # genes.bed name col = "gene_name|gene_id|gene_type"
    _biotype = {}
    for cols in read_bed(genes_bed):
        if len(cols) < 4:
            continue
        parts = cols[3].split("|")
        if len(parts) >= 3:
            _biotype[parts[0]] = parts[2]

    fusions = []

    # ── Optionally expand breakpoint intervals for gene intersection ──────────
    if bp_window > 0:
        _exp = _tf.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
        try:
            for cols in read_bed(bp_bed):
                chrom, start, end, bp_id, svtype = cols[:5]
                s = max(0, int(start) - bp_window)
                e = int(end) + bp_window
                _exp.write(f"{chrom}\t{s}\t{e}\t{bp_id}\t{svtype}\n")
            _exp.close()
            _bp = _exp.name
        except Exception:
            os.unlink(_exp.name)
            raise
    else:
        _bp = bp_bed

    # ── Intrachromosomal ────────────────────────────────────────────────────
    # Collect ALL genes overlapping each (possibly expanded) breakpoint.
    # When genes are nested (e.g. CHIC2 inside KIT), a single-gene selector
    # may drop the clinically relevant partner. All left × right combinations
    # are checked so both CHIC2::PDGFRA and KIT::PDGFRA are considered.

    # bp_id → [(gene_name, element), ...]
    bp_all = defaultdict(list)

    # exon hits: group by (bp_id, gene_name) to resolve best exon element
    _exon = defaultdict(list)   # bp_id → [(gene, exon_num, gtype)]
    for cols in intersect(_bp, exons_detail_bed):
        if len(cols) < 9:
            continue
        bp_id = cols[3]
        parts = cols[8].split("|")
        gene_name = parts[0]
        exon_num  = parts[2] if len(parts) > 2 else "?"
        gtype     = parts[3] if len(parts) > 3 else "unknown"
        _exon[bp_id].append((gene_name, exon_num, gtype))

    # gene body hits: all genes per bp_id
    _gbody = defaultdict(set)   # bp_id → {gene_name}
    for cols in intersect(_bp, genes_bed):
        if len(cols) < 9:
            continue
        _gbody[cols[3]].add(gene_name_from_field(cols[8]))

    if bp_window > 0:
        os.unlink(_bp)

    # For each gene overlapping a BP, resolve its element
    for bp_id in set(list(_exon.keys()) + list(_gbody.keys())):
        exon_by_gene = defaultdict(list)
        for gname, enum, gtype in _exon.get(bp_id, []):
            exon_by_gene[gname].append((enum, gtype))

        for gname, hits in exon_by_gene.items():
            try:
                best_enum = min(hits, key=lambda x: int(x[0]))[0]
            except ValueError:
                best_enum = hits[0][0]
            bp_all[bp_id].append((gname, f"exon_{best_enum}"))

        for gname in _gbody.get(bp_id, set()):
            if gname not in exon_by_gene:
                bp_all[bp_id].append((gname, "intronic"))

    # Build per-SV breakpoint data
    sv_bp = {}
    for cols in read_bed(bp_bed):
        chrom, start, _, bp_id, svtype = cols[:5]
        sv_id, side = bp_id.rsplit("|", 1)
        sv_bp.setdefault(sv_id, {"chrom": chrom, "svtype": svtype})
        sv_bp[sv_id][side] = {"genes": bp_all.get(bp_id, []), "pos": start}

    seen_intra = set()
    for sv_id, data in sv_bp.items():
        left_genes  = data.get("left",  {}).get("genes", [])
        right_genes = data.get("right", {}).get("genes", [])
        left_pos    = data.get("left",  {}).get("pos", "")
        right_pos   = data.get("right", {}).get("pos", "")
        if not left_genes or not right_genes:
            continue

        # Skip if both breakpoints are within the same overlapping region
        # (tiny SV — cannot span between two distinct gene bodies)
        try:
            sv_size = abs(int(right_pos) - int(left_pos))
        except (ValueError, TypeError):
            sv_size = 0
        if sv_size < min_fusion_size:
            continue

        for ga, ea in left_genes:
            for gb, eb in right_genes:
                if not ga or not gb or ga == gb:
                    continue
                pair_key = (sv_id, tuple(sorted([ga, gb])))
                if pair_key in seen_intra:
                    continue
                seen_intra.add(pair_key)
                any_in  = gene_set and (ga in gene_set or gb in gene_set)
                both_in = gene_set and (ga in gene_set and gb in gene_set)
                ta = _biotype.get(ga, "unknown")
                tb = _biotype.get(gb, "unknown")
                support, af = sv_info.get(sv_id, (".", "."))
                fusions.append([
                    sv_id, f"{ga}::{gb}", data["svtype"], "intrachromosomal",
                    data["chrom"], left_pos,  ga, ea, ta,
                    data["chrom"], right_pos, gb, eb, tb,
                    "yes" if ea.startswith("exon_") and eb.startswith("exon_") else "no",
                    "yes" if ta == "protein_coding" and tb == "protein_coding" else "no",
                    "yes" if any_in else "no",
                    "yes" if both_in else "no",
                    support, af,
                ])

    # ── BND translocations ───────────────────────────────────────────────────
    # Use own_annot (BND's position) and mate_annot (ALT-encoded mate position).
    # This works even when only one BND record per pair is in the VCF.
    # Deduplicate with a set of sorted gene pairs to avoid double-counting
    # when both mates ARE present.
    own_annot, mate_annot, bnd_meta = bnd_context
    seen_pairs = set()

    for bnd_id, rec in bnd_meta.items():
        if rec["mate_chrom"] in (None, ".") or rec["mate_pos"] is None:
            continue
        ga, ea = own_annot.get(bnd_id,  ("", "intergenic"))
        gb, eb = mate_annot.get(bnd_id, ("", "intergenic"))
        if not ga or not gb or ga == gb:
            continue
        pair_key = tuple(sorted([ga, gb]))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)
        any_in  = gene_set and (ga in gene_set or gb in gene_set)
        both_in = gene_set and (ga in gene_set and gb in gene_set)
        ta = _biotype.get(ga, "unknown")
        tb = _biotype.get(gb, "unknown")
        support, af = sv_info.get(bnd_id, (".", "."))
        fusions.append([
            bnd_id, f"{ga}::{gb}", "BND", "translocation",
            rec["chrom"], str(rec["pos"] - 1), ga, ea, ta,
            rec["mate_chrom"], str(rec["mate_pos"] - 1), gb, eb, tb,
            "yes" if ea.startswith("exon_") and eb.startswith("exon_") else "no",
            "yes" if ta == "protein_coding" and tb == "protein_coding" else "no",
            "yes" if any_in else "no",
            "yes" if both_in else "no",
            support, af,
        ])

    header = (
        "sv_id\tfusion_name\tsvtype\tsource\t"
        "chrom_A\tpos_A\tgene_A\telement_A\tgene_type_A\t"
        "chrom_B\tpos_B\tgene_B\telement_B\tgene_type_B\t"
        "both_exonic\tboth_protein_coding\t"
        "any_in_gene_list\tboth_in_gene_list\t"
        "support\taf\n"
    )
    with open(outfile, "w") as f:
        f.write(header)
        for row in fusions:
            f.write("\t".join(str(x) for x in row) + "\n")
    return len(fusions)


def load_gene_list(path):
    """Return a set of gene names from a file (handles optional leading line numbers)."""
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # accept both "GENE" and "123  GENE" formats
            gene = parts[1] if len(parts) >= 2 and parts[0].isdigit() else parts[0]
            genes.add(gene)
    return genes


def write_fusions_gbm_protein(fusions_tsv, outfile):
    """Write fusions where any_in_gene_list=yes AND both genes are protein_coding."""
    hits = 0
    with open(fusions_tsv) as src, open(outfile, "w") as out:
        header = next(src)
        cols_hdr = header.rstrip("\n").split("\t")
        try:
            idx_any = cols_hdr.index("any_in_gene_list")
            idx_ta  = cols_hdr.index("gene_type_A")
            idx_tb  = cols_hdr.index("gene_type_B")
        except ValueError as e:
            sys.exit(f"Column not found in {fusions_tsv}: {e}")
        out.write(header)
        for line in src:
            cols = line.rstrip("\n").split("\t")
            if (len(cols) > max(idx_any, idx_ta, idx_tb)
                    and cols[idx_any] == "yes"
                    and cols[idx_ta] == "protein_coding"
                    and cols[idx_tb] == "protein_coding"):
                out.write(line)
                hits += 1
    return hits


def write_fusions_filtered(fusions_tsv, col, outfile):
    """Write fusions.tsv rows where column `col` == 'yes'."""
    hits = 0
    with open(fusions_tsv) as src, open(outfile, "w") as out:
        header = next(src)
        cols_hdr = header.rstrip("\n").split("\t")
        try:
            idx = cols_hdr.index(col)
        except ValueError:
            sys.exit(f"Column '{col}' not found in {fusions_tsv}")
        out.write(header)
        for line in src:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > idx and cols[idx] == "yes":
                out.write(line)
                hits += 1
    return hits


def build_bnd_context(bnd_bed, exons_detail_bed, genes_bed, vcf, window=0):
    """Annotate BND positions and mate positions; parse ALT field once.

    Returns (own_annot, mate_annot, bnd_meta) shared by annotate_bnd and
    detect_fusions.

    own_annot  : {bnd_id: (gene, element)} — annotation of the BND's own position
    mate_annot : {bnd_id: (gene, element)} — annotation of the ALT-encoded mate pos
    bnd_meta   : {bnd_id: {chrom, pos, mate_chrom, mate_pos}}

    Both own and mate positions are expanded by `window` bp before intersecting
    with gene/exon annotations.  This approach does NOT require the mate BND
    record to be present in the VCF — critical for wf-human-variation output
    which keeps only one record per BND pair.
    """
    import tempfile as _tf

    # ── Parse ALT field to get mate positions ────────────────────────────────
    raw = run(
        f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%ALT\\n'"
        f" -i 'SVTYPE=\"BND\"' {vcf}",
        check=False,
    )
    bnd_meta = {}
    for line in raw.splitlines():
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        chrom, pos, bnd_id, alt = cols
        try:
            pos_i = int(pos)
        except ValueError:
            continue
        mc, mp = parse_bnd_alt(alt)
        bnd_meta[bnd_id] = {"chrom": chrom, "pos": pos_i,
                             "mate_chrom": mc or ".", "mate_pos": mp}

    # ── Annotate own positions (from bnd_bed, already 1-bp) ─────────────────
    own_exp = _tf.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    try:
        for cols in read_bed(bnd_bed):
            chrom, start, end, rec_id, svtype = cols[:5]
            s = max(0, int(start) - window)
            e = int(end) + window
            own_exp.write(f"{chrom}\t{s}\t{e}\t{rec_id}\t{svtype}\n")
        own_exp.close()
        own_annot = annotate_positions(own_exp.name, exons_detail_bed, genes_bed)
    finally:
        os.unlink(own_exp.name)

    # ── Annotate mate positions (from ALT field) ─────────────────────────────
    mate_exp = _tf.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    try:
        for bnd_id, rec in bnd_meta.items():
            mc, mp = rec["mate_chrom"], rec["mate_pos"]
            if mc == "." or mp is None:
                continue
            s = max(0, mp - 1 - window)  # mp is 1-based VCF POS
            e = mp + window
            mate_exp.write(f"{mc}\t{s}\t{e}\t{bnd_id}\t.\n")
        mate_exp.close()
        mate_annot = annotate_positions(mate_exp.name, exons_detail_bed, genes_bed)
    finally:
        os.unlink(mate_exp.name)

    return own_annot, mate_annot, bnd_meta


def annotate_bnd(bnd_bed, bnd_context, outfile):
    """Write sv_bnd.tsv using pre-built BND context.

    Columns: sv_id, chrom, pos, svtype, gene, element,
             mate_chrom, mate_pos, mate_gene, mate_element, fusion_name

    Works even when the mate BND record is absent from the VCF.
    """
    own_annot, mate_annot, bnd_meta = bnd_context

    with open(bnd_bed) as bf, open(outfile, "w") as out:
        out.write("sv_id\tchrom\tpos\tsvtype\tgene\telement\t"
                  "mate_chrom\tmate_pos\tmate_gene\tmate_element\tfusion_name\n")
        for line in bf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, start, _, sv_id, svtype = cols[:5]
            gene, element = own_annot.get(sv_id, ("", "intergenic"))
            rec = bnd_meta.get(sv_id, {})
            mc  = rec.get("mate_chrom", ".")
            mp  = rec.get("mate_pos")
            mg, me = mate_annot.get(sv_id, ("", "intergenic"))
            fname = f"{gene}::{mg}" if gene and mg and gene != mg else ""
            out.write(f"{sv_id}\t{chrom}\t{start}\t{svtype}\t"
                      f"{gene}\t{element}\t"
                      f"{mc}\t{mp or '.'}\t{mg}\t{me}\t{fname}\n")


def write_breakpoints_filtered(bp_tsv, gene_set, outfile):
    """Filter sv_breakpoints.tsv to rows whose gene is in gene_set."""
    hits = 0
    with open(bp_tsv) as src, open(outfile, "w") as out:
        header = next(src)
        out.write(header)
        for line in src:
            cols = line.rstrip("\n").split("\t")
            # gene is col index 5 (sv_id, bp_side, chrom, bp_pos, svtype, gene, element)
            if len(cols) >= 6 and cols[5] in gene_set:
                out.write(line)
                hits += 1
    return hits


def write_gene_filtered(results, gene_set, outfile):
    """Write SVs that overlap at least one gene in gene_set, adding 'gbm_genes' column."""
    hits = 0
    with open(outfile, "w") as out:
        out.write("sv_id\tchrom\tstart\tend\tsvtype\tclassification\tgenes\tgbm_genes\n")
        for row in results:
            sv_id, chrom, start, end, svtype, cls, genes = row
            if not genes:
                continue
            sv_genes = set(genes.split(","))
            matched = sorted(sv_genes & gene_set)
            if matched:
                out.write("\t".join(row) + "\t" + ",".join(matched) + "\n")
                hits += 1
    return hits


def print_summary(classified):
    from collections import Counter
    counts = Counter(r[5] for r in classified)
    print("\n=== Classification summary ===", file=sys.stderr)
    for cls in ("full_gene_loss", "exon_disruption", "intronic_disruption",
                "intergenic", "unclassified_bnd"):
        print(f"  {cls:25s}: {counts.get(cls, 0):>6,}", file=sys.stderr)
    print(f"  {'TOTAL':25s}: {sum(counts.values()):>6,}", file=sys.stderr)


def sample_id_from_vcf(vcf_path):
    """Derive a sample ID from the VCF filename (strips .vcf.gz / .vcf)."""
    base = os.path.basename(vcf_path)
    for ext in (".vcf.gz", ".vcf.bgz", ".vcf"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base.split(".")[0]


def parse_sample_list(path):
    """Read sample IDs from a file (one per line, # lines are comments)."""
    samples = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                samples.append(line)
    return samples


def process_sample(vcf, sid, outdir, args, genes_bed, exons_bed,
                   exons_detail_bed, gene_set):
    """Run the full pipeline for one sample using pre-built annotation BEDs."""

    def outpath(suffix):
        return os.path.join(outdir, f"{sid}.{suffix}")

    out_sv              = outpath("sv_classified.tsv")
    out_bnd             = outpath("sv_bnd.tsv")
    out_bp              = outpath("sv_breakpoints.tsv")
    out_fusions         = outpath("sv_fusions.tsv")
    out_fusions_any_any = outpath("sv_fusions_any_any.tsv")
    out_fusions_any     = outpath("sv_fusions_any_gbm.tsv")
    out_fusions_both    = outpath("sv_fusions_both_gbm.tsv")
    out_fusions_protein = outpath("sv_fusions_gbm_protein.tsv")
    out_genes           = outpath("sv_gbm_genes.tsv")
    out_bp_genes        = outpath("sv_breakpoints_gbm.tsv")

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Sample : {sid}", file=sys.stderr)
    print(f"VCF    : {vcf}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmp:
        sv_bed  = os.path.join(tmp, "sv.bed")
        bnd_bed = os.path.join(tmp, "bnd.bed")

        print("Extracting SVs from VCF...", file=sys.stderr)
        extract_sv_bed(vcf, sv_bed, bnd_bed)

        print("Running bedtools intersections...", file=sys.stderr)
        fgl_lines = intersect(genes_bed, sv_bed, "-f 1.0")
        full_gene_loss  = build_map_gene_a(fgl_lines, sv_id_col=9, gene_col=3)
        exon_lines      = intersect(sv_bed, exons_bed)
        exon_disruption = build_map_sv_a(exon_lines, sv_id_col=3, gene_col=8)
        gene_lines      = intersect(sv_bed, genes_bed)
        gene_overlap    = build_map_sv_a(gene_lines, sv_id_col=3, gene_col=8)

        print("Classifying SVs...", file=sys.stderr)
        results = []
        for cols in read_bed(sv_bed):
            chrom, start, end, sv_id, svtype = cols[:5]
            if sv_id in full_gene_loss:
                cls, genes = "full_gene_loss", ",".join(sorted(full_gene_loss[sv_id]))
            elif sv_id in exon_disruption:
                cls, genes = "exon_disruption", ",".join(sorted(exon_disruption[sv_id]))
            elif sv_id in gene_overlap:
                cls, genes = "intronic_disruption", ",".join(sorted(gene_overlap[sv_id]))
            else:
                cls, genes = "intergenic", ""
            results.append([sv_id, chrom, start, end, svtype, cls, genes])

        for cols in read_bed(bnd_bed):
            chrom, start, end, sv_id, svtype = cols[:5]
            results.append([sv_id, chrom, start, end, svtype, "unclassified_bnd", ""])

        with open(out_sv, "w") as fh:
            fh.write("sv_id\tchrom\tstart\tend\tsvtype\tclassification\tgenes\n")
            for row in results:
                fh.write("\t".join(row) + "\n")

        bp_bed = os.path.join(tmp, "breakpoints.bed")
        make_breakpoints_bed(sv_bed, bp_bed)
        print("Annotating breakpoints...", file=sys.stderr)
        annotate_breakpoints(bp_bed, exons_detail_bed, genes_bed, out_bp)

        print("Annotating BND breakpoints...", file=sys.stderr)
        bnd_context = build_bnd_context(bnd_bed, exons_detail_bed, genes_bed,
                                        vcf, window=args.bnd_window)
        annotate_bnd(bnd_bed, bnd_context, out_bnd)

        print("Detecting gene fusions...", file=sys.stderr)
        sv_info = extract_sv_info(vcf)
        n = detect_fusions(bp_bed, bnd_context, exons_detail_bed, genes_bed,
                           gene_set, out_fusions,
                           min_fusion_size=args.min_fusion_size,
                           bp_window=args.bp_window,
                           sv_info=sv_info)
        print(f"Fusion candidates: {n}", file=sys.stderr)

    print_summary(results)

    # any::any — all fusions regardless of gene list
    shutil.copy(out_fusions, out_fusions_any_any)

    if gene_set:
        hits = write_gene_filtered(results, gene_set, out_genes)
        print(f"GBM-filtered SVs: {hits}", file=sys.stderr)
        write_breakpoints_filtered(out_bp, gene_set, out_bp_genes)
        write_fusions_filtered(out_fusions, "any_in_gene_list",  out_fusions_any)
        write_fusions_filtered(out_fusions, "both_in_gene_list", out_fusions_both)
        n = write_fusions_gbm_protein(out_fusions, out_fusions_protein)
        print(f"GBM + protein_coding fusions: {out_fusions_protein} ({n})",
              file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Classify Sniffles SVs by gene impact — single or multi-sample",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # single sample
  python classify_sv.py --vcf T25-000.vcf.gz --gff3 gencode.v48.annotation.gff3 \\
      --gene-list occ_fusions_genes.txt --outdir results/

  # multiple samples from a list
  python classify_sv.py --sample-list samples.txt --input-dir /data/vcf/ \\
      --gff3 gencode.v48.annotation.gff3 --gene-list occ_fusions_genes.txt \\
      --outdir results/

samples.txt — one sample ID per line (# = comment):
  T25-000
  # 
""")

    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--vcf", default=None,
                      help="Single-sample VCF (bgzipped+tabixed or plain)")
    mode.add_argument("--sample-list", default=None,
                      help="File with one sample ID per line")

    ap.add_argument("--gff3", required=True, help="GENCODE GFF3 annotation")
    ap.add_argument("--input-dir", default=None,
                    help="Directory containing VCF files (multi-sample mode)")
    ap.add_argument("--vcf-suffix", default=".sniffles.vcf.gz",
                    help="Suffix appended to each sample ID to form its VCF filename "
                         "(default: .sniffles.vcf.gz). "
                         "VCF path = {input_dir}/{sample_id}{vcf_suffix}")
    ap.add_argument("--sample-id", default=None,
                    help="Override sample ID (single-sample mode only; "
                         "default: derived from VCF filename)")
    ap.add_argument("--outdir", default=".",
                    help="Output directory (default: current directory)")
    ap.add_argument("--gene-list", default=None,
                    help="Gene list file (e.g. GBM drivers); one gene per line, "
                         "optional leading line number")
    ap.add_argument("--protein-coding-only", action="store_true",
                    help="Restrict gene/exon intersections to protein_coding genes")
    ap.add_argument("--bnd-window", type=int, default=1000,
                    help="Window (bp) around BND positions for gene annotation "
                         "(default: 1000)")
    ap.add_argument("--bp-window", type=int, default=0,
                    help="Window (bp) added around each intrachromosomal breakpoint "
                         "before gene intersection in fusion detection (default: 0). "
                         "Increase (e.g. 200) to recover fusions where Sniffles "
                         "places a breakpoint just outside a gene boundary. "
                         "Does NOT affect sv_breakpoints.tsv — only fusion detection.")
    ap.add_argument("--min-fusion-size", type=int, default=10_000,
                    help="Minimum SV size (bp) for intrachromosomal fusions "
                         "(default: 10000). Filters tiny DELs where both breakpoints "
                         "lie within overlapping nested genes — these are annotation "
                         "artifacts, not real fusions.")
    args = ap.parse_args()

    # ── Validate arguments ────────────────────────────────────────────────────
    if args.sample_list and not args.input_dir:
        ap.error("--input-dir is required when using --sample-list")

    os.makedirs(args.outdir, exist_ok=True)
    gene_types = {"protein_coding"} if args.protein_coding_only else None
    gene_set   = load_gene_list(args.gene_list) if args.gene_list else None

    # ── Collect samples ───────────────────────────────────────────────────────
    if args.vcf:
        sid = args.sample_id or sample_id_from_vcf(args.vcf)
        samples = [(sid, args.vcf)]
    else:
        sample_ids = parse_sample_list(args.sample_list)
        samples = []
        for sid in sample_ids:
            vcf = os.path.join(args.input_dir, f"{sid}{args.vcf_suffix}")
            if not os.path.exists(vcf):
                print(f"Warning: VCF not found for {sid}: {vcf} — skipping",
                      file=sys.stderr)
                continue
            samples.append((sid, vcf))
        print(f"Found {len(samples)}/{len(sample_ids)} VCFs from "
              f"{args.sample_list}", file=sys.stderr)

    # ── Build annotation BEDs once — shared across all samples ───────────────
    ann_tmp = tempfile.mkdtemp(prefix="classify_sv_ann_")
    try:
        genes_bed        = os.path.join(ann_tmp, "genes.bed")
        exons_bed        = os.path.join(ann_tmp, "exons.bed")
        exons_detail_bed = os.path.join(ann_tmp, "exons_detail.bed")
        print("Parsing GENCODE GFF3 (once for all samples)...", file=sys.stderr)
        extract_annotation_beds(args.gff3, genes_bed, exons_bed, exons_detail_bed,
                                gene_types)

        # ── Process each sample ───────────────────────────────────────────────
        for i, (sid, vcf) in enumerate(samples, 1):
            print(f"\n[{i}/{len(samples)}]", file=sys.stderr)
            process_sample(vcf, sid, args.outdir, args,
                           genes_bed, exons_bed, exons_detail_bed, gene_set)

    finally:
        shutil.rmtree(ann_tmp, ignore_errors=True)

    print("\nAll done.", file=sys.stderr)


if __name__ == "__main__":
    main()
