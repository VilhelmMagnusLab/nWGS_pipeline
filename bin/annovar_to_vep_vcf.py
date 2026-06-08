#!/usr/bin/env python3
"""
Convert ANNOVAR-format ClinVar/COSMIC annotation files to bgzipped,
tabix-indexed VCF for use as VEP --custom files.

Usage:
  annovar_to_vep_vcf.py clinvar <hg38_clinvar.txt> <output.vcf.gz>
  annovar_to_vep_vcf.py cosmic  <hg38_cosmic.txt>  <output.vcf.gz>

ANNOVAR format: #Chr Start End Ref Alt [fields...]
VEP --custom expects: bgzipped VCF, tabix-indexed, chr-prefixed chromosomes.
"""

import sys
import subprocess
import os
import tempfile


CHROM_ORDER = (
    [f"chr{i}" for i in range(1, 23)]
    + ["chrX", "chrY", "chrM", "chrMT"]
)


def sanitize(val):
    return val.replace(";", ",").replace(" ", "_").replace("\t", "_") if val else "."


def convert_clinvar(txt_path, vcf_gz_path):
    header_lines = [
        "##fileformat=VCFv4.1",
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance">',
        '##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease name">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    # Columns: Chr Start End Ref Alt CLNALLELEID CLNDN CLNDISDB CLNREVSTAT CLNSIG
    records = []
    with open(txt_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            chrom, pos, ref, alt = c[0], c[1], c[3], c[4]
            if ref == "-" or alt == "-":
                continue
            chrom = chrom if chrom.startswith("chr") else "chr" + chrom
            clnsig = sanitize(c[9])
            clndn = sanitize(c[6])
            info = f"CLNSIG={clnsig};CLNDN={clndn}"
            records.append((chrom, int(pos), ref, alt, info))

    _write_sorted_vcf(header_lines, records, vcf_gz_path)


def convert_cosmic(txt_path, vcf_gz_path):
    header_lines = [
        "##fileformat=VCFv4.1",
        '##INFO=<ID=GENE,Number=.,Type=String,Description="COSMIC gene / variant ID">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    # Columns: Chr Start End Ref Alt COSMIC100  (value like "ID=COSV70831266")
    records = []
    with open(txt_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            chrom, pos, ref, alt = c[0], c[1], c[3], c[4]
            if ref == "-" or alt == "-":
                continue
            chrom = chrom if chrom.startswith("chr") else "chr" + chrom
            # Strip "ID=" prefix so INFO becomes GENE=COSV...
            cosmic_val = sanitize(c[5]).replace("ID=", "")
            info = f"GENE={cosmic_val}" if cosmic_val and cosmic_val != "." else "GENE=."
            records.append((chrom, int(pos), ref, alt, info))

    _write_sorted_vcf(header_lines, records, vcf_gz_path)


def _write_sorted_vcf(header_lines, records, vcf_gz_path):
    chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
    records.sort(key=lambda r: (chrom_rank.get(r[0], 999), r[1]))

    tmp_vcf = vcf_gz_path.replace(".vcf.gz", "_tmp.vcf")
    with open(tmp_vcf, "w") as out:
        out.write("\n".join(header_lines) + "\n")
        for chrom, pos, ref, alt, info in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    subprocess.run(["bgzip", "-f", tmp_vcf], check=True)
    bgz = tmp_vcf + ".gz"
    os.rename(bgz, vcf_gz_path)
    subprocess.run(["tabix", "-p", "vcf", vcf_gz_path], check=True)
    print(f"Written: {vcf_gz_path} (+.tbi)")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("clinvar", "cosmic"):
        sys.exit("Usage: annovar_to_vep_vcf.py clinvar|cosmic <input.txt> <output.vcf.gz>")
    mode, inp, out = sys.argv[1], sys.argv[2], sys.argv[3]
    if mode == "clinvar":
        convert_clinvar(inp, out)
    else:
        convert_cosmic(inp, out)
