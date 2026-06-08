#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a VEP-annotated VCF to ANNOVAR-equivalent TSV for the DIANA pipeline.

Produces the column layout expected by merge_annotations_prospective.R:

  clair3   mode — 19 cols; Otherinfo13 (col 19) = SAMPLE
  clairsto mode — 18 cols; Otherinfo10 (col 18) = SAMPLE

Usage:
  vep_to_annovar_format.py <vep_annotated.vcf[.gz]> <output.tsv> clair3|clairsto
"""

from __future__ import print_function
import sys
import re
import gzip
import io

# ── consequence → ANNOVAR Func.refGene ──────────────────────────────────────
FUNC = {
    'upstream_gene_variant':              'upstream',
    'downstream_gene_variant':            'downstream',
    'missense_variant':                   'exonic',
    'synonymous_variant':                 'exonic',
    'stop_gained':                        'exonic',
    'stop_lost':                          'exonic',
    'frameshift_variant':                 'exonic',
    'inframe_insertion':                  'exonic',
    'inframe_deletion':                   'exonic',
    'start_lost':                         'exonic',
    'protein_altering_variant':           'exonic',
    'splice_donor_variant':               'splicing',
    'splice_acceptor_variant':            'splicing',
    'splice_region_variant':              'splicing',
    'intron_variant':                     'intronic',
    'intergenic_variant':                 'intergenic',
    'non_coding_transcript_exon_variant': 'ncRNA_exonic',
    'non_coding_transcript_variant':      'ncRNA_intronic',
    '3_prime_UTR_variant':                'UTR3',
    '5_prime_UTR_variant':                'UTR5',
}

# ── consequence → ANNOVAR ExonicFunc.refGene ─────────────────────────────────
EXONIC = {
    'missense_variant':   'nonsynonymous SNV',
    'synonymous_variant': 'synonymous SNV',
    'stop_gained':        'stopgain',
    'stop_lost':          'stoploss',
    'frameshift_variant': 'frameshift substitution',
    'inframe_insertion':  'nonframeshift insertion',
    'inframe_deletion':   'nonframeshift deletion',
    'start_lost':         'nonsynonymous SNV',
}

# Severity order used to pick one consequence when multiple are present
PRIORITY = [
    'frameshift_variant', 'stop_gained', 'start_lost', 'stop_lost',
    'splice_donor_variant', 'splice_acceptor_variant',
    'inframe_insertion', 'inframe_deletion', 'missense_variant',
    'splice_region_variant', 'synonymous_variant',
    'upstream_gene_variant', 'downstream_gene_variant',
    'intron_variant', '3_prime_UTR_variant', '5_prime_UTR_variant',
    'non_coding_transcript_exon_variant', 'non_coding_transcript_variant',
    'intergenic_variant',
]


def open_file(path):
    if path.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'))
    return open(path)


def primary_consequence(csq_str):
    terms = csq_str.split('&')
    for p in PRIORITY:
        if p in terms:
            return p
    return terms[0] if terms else ''


def parse_csq_format(header_line):
    m = re.search(r'Format: ([^"]+)"', header_line)
    return m.group(1).split('|') if m else []


def compute_end(pos, ref, alt):
    pos = int(pos)
    return str(pos + len(ref) - 1) if len(ref) >= len(alt) else str(pos)


def build_aa_change(csq):
    """Reconstruct ANNOVAR-style GENE:TRANSCRIPT:exonN:HGVSc:HGVSp string."""
    symbol  = csq.get('SYMBOL', '')
    feature = csq.get('Feature', '').split('.')[0]  # strip version (NM_002167.5 → NM_002167)
    exon    = csq.get('EXON', '')
    hgvsc   = csq.get('HGVSc', '')
    hgvsp   = csq.get('HGVSp', '')

    parts = [symbol, feature]
    if exon:
        parts.append('exon' + exon.split('/')[0])
    if hgvsc:
        parts.append(hgvsc.split(':')[-1] if ':' in hgvsc else hgvsc)
    if hgvsp:
        parts.append(hgvsp.split(':')[-1] if ':' in hgvsp else hgvsp)
    parts = [p for p in parts if p]
    return ':'.join(parts) if len(parts) > 1 else ''


def convert(vcf_path, out_path, mode='clair3'):
    if mode == 'clair3':
        # 19 cols — matches ANNOVAR cut -f1-16,26,28,29 from 29-col multianno
        header = [
            'Chr', 'Start', 'End', 'Ref', 'Alt',
            'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
            'ExonicFunc.refGene', 'AAChange.refGene',
            'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG',
            'COSMIC100',
            'Otherinfo10',   # FILTER  (col 17)
            'Otherinfo12',   # FORMAT  (col 18)
            'Otherinfo13',   # SAMPLE  (col 19) — parsed by R parse_format_field
        ]
    else:
        # 18 cols — matches expected layout for clairsto (consistent header+data)
        header = [
            'Chr', 'Start', 'End', 'Ref', 'Alt',
            'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
            'ExonicFunc.refGene', 'AAChange.refGene',
            'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG',
            'Otherinfo1',    # FILTER  (col 16)
            'Otherinfo9',    # FORMAT  (col 17)
            'Otherinfo10',   # SAMPLE  (col 18) — parsed by R parse_format_field
        ]

    csq_fields = []

    with open_file(vcf_path) as fin, open(out_path, 'w') as fout:
        fout.write('\t'.join(header) + '\n')

        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('##'):
                if 'ID=CSQ' in line:
                    csq_fields = parse_csq_format(line)
                continue
            if line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) < 8:
                continue

            chrom  = cols[0]
            pos    = cols[1]
            ref    = cols[3].upper()
            alt    = cols[4].split(',')[0].upper()
            filt   = cols[6]
            info   = cols[7]
            fmt    = cols[8] if len(cols) > 8 else ''
            sample = cols[9] if len(cols) > 9 else ''

            end = compute_end(pos, ref, alt)

            # Parse first CSQ entry (--pick gives one per variant)
            csq = {}
            m = re.search(r'CSQ=([^;]+)', info)
            if m and csq_fields:
                vals = m.group(1).split(',')[0].split('|')
                for i, name in enumerate(csq_fields):
                    csq[name] = vals[i] if i < len(vals) else ''

            pcsq        = primary_consequence(csq.get('Consequence', ''))
            func        = FUNC.get(pcsq, 'intronic')
            gene        = csq.get('SYMBOL') or csq.get('Gene', '')
            distance    = csq.get('DISTANCE', '')
            gene_detail = ('dist=' + distance) if (distance and func == 'upstream') else ''
            exonic_func = EXONIC.get(pcsq, '') if func == 'exonic' else ''
            aa_change   = build_aa_change(csq)
            clnsig      = csq.get('ClinVar_CLNSIG', '')
            clndn       = csq.get('ClinVar_CLNDN', '')
            clndisdb    = csq.get('ClinVar_CLNDISDB', '')
            clnrevstat  = csq.get('ClinVar_CLNREVSTAT', '')
            cosmic_raw  = csq.get('COSMIC_GENE', '')
            cosmic100   = ('ID=' + cosmic_raw) if cosmic_raw else ''

            if mode == 'clair3':
                row = [
                    chrom, pos, end, ref, alt,
                    func, gene, gene_detail, exonic_func, aa_change,
                    '', clndn, clndisdb, clnrevstat, clnsig,
                    cosmic100,
                    filt, fmt, sample,
                ]
            else:
                row = [
                    chrom, pos, end, ref, alt,
                    func, gene, gene_detail, exonic_func, aa_change,
                    '', clndn, clndisdb, clnrevstat, clnsig,
                    filt, fmt, sample,
                ]

            fout.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('Usage: ' + sys.argv[0] + ' <vep.vcf[.gz]> <output.tsv> clair3|clairsto')
    convert(sys.argv[1], sys.argv[2], sys.argv[3])
