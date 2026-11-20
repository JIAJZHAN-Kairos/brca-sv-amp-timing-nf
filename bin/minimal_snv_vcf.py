#!/usr/bin/env python3
"""Rewrite a SASH somatic VCF as a minimal VCF that SVclone's pyvcf 'mutect' loader can parse:
PASS biallelic SNVs only, FORMAT reduced to GT:AD for the tumour sample (column named <sample>).
Some SASH VCFs carry FORMAT fields declared Number=1 with several values, which pyvcf rejects.
Usage: minimal_snv_vcf.py <in.vcf.gz> <sample> <out.vcf>"""
import gzip, sys
src, sample, dst = sys.argv[1:4]
n = 0
with gzip.open(src, "rt") as f, open(dst, "w") as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
    for ln in f:
        if ln.startswith("##contig"):
            o.write(ln)
        if ln.startswith("##"):
            continue
        c = ln.rstrip("\n").split("\t")
        if ln.startswith("#CHROM"):
            ti = c.index(sample)
            o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
            continue
        if c[6] != "PASS" or len(c[3]) != 1 or len(c[4]) != 1:
            continue
        fmt = c[8].split(":"); val = c[ti].split(":")
        ad = val[fmt.index("AD")].split(",")
        if len(ad) != 2:
            continue
        o.write("\t".join(c[:5] + [".", "PASS", ".", "GT:AD", "0/1:%s,%s" % (ad[0], ad[1])]) + "\n")
        n += 1
print(f"minimal_snv_vcf: {n} PASS SNVs for {sample}", file=sys.stderr)
