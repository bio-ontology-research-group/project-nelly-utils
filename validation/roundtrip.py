#!/usr/bin/env python3
"""End-to-end simulator validation for ONE disease, on a single chromosome.

  simulate (haplotype round-trip) -> align to hg38 -> call -> compare to truth.

Modes:
  buggy  -- uses scripts/map_and_simulate.py's linear-interpolation target pos
  fixed  -- uses the indel-aware CIGAR-walk target pos

Reports whether the planted variant is recovered at the EXACT hg38 position and
whether it is called HOMOZYGOUS (recessive) -- the two things that were wrong.
"""
import os, sys, json, subprocess, tempfile
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confirm_drift import (get_context, map_paf, buggy_target,
                           cigar_walk_target, HAP1, HAP2)
# 'fixed' mode validates the ACTUAL shipped fix, not a copy. Stub the optional
# phenotype/phenopackets dependency (not needed for coordinate logic).
import types
_ph = types.ModuleType("phenotype")
_ph.create_phenopacket = lambda *a, **k: None
sys.modules.setdefault("phenotype", _ph)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from map_and_simulate import calculate_target_pos as src_calculate_target_pos

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(BASE, "config", "diseases.json")
WORK = os.path.join(BASE, "validation", "work")
ROI = 2500          # +/- bp extracted per haplotype for read simulation
COV = 30            # per-haplotype coverage (art); diploid ~= 2x at the locus


def rc(s):
    c = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
         'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(c.get(b, b) for b in reversed(s))


def extract_patch(fa_path, contig, target_pos, ref, alt, strand, do_patch):
    with pysam.FastaFile(fa_path) as fa:
        clen = fa.get_reference_length(contig)
        start = max(0, target_pos - 1 - ROI)
        end = min(clen, target_pos - 1 + len(ref) + ROI)
        seq = fa.fetch(contig, start, end)
    rel = (target_pos - 1) - start
    p_ref, p_alt = (rc(ref), rc(alt)) if strand == '-' else (ref, alt)
    observed = seq[rel:rel + len(p_ref)].upper()
    status = "OK" if observed == p_ref.upper() else f"REFMISMATCH(saw {observed} want {p_ref})"
    if do_patch:
        seq = seq[:rel] + p_alt + seq[rel + len(p_ref):]
    return seq, status


def run(disease, mode):
    diseases = json.load(open(CONFIG))
    d = next(x for x in diseases if disease.lower() in x["name"].lower())
    v = d["variant"]
    chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]
    recessive = d.get("inheritance") == "recessive"
    print(f"\n=== {d['name']}  {chrom}:{pos} {ref}>{alt}  "
          f"({'recessive->HOM' if recessive else 'dominant->HET'})  mode={mode} ===")

    seq, voff = get_context(chrom, pos, ref)
    sample_dir = os.path.join(WORK, f"{disease.lower()}_{mode}")
    os.makedirs(sample_dir, exist_ok=True)
    roi_fa = os.path.join(sample_dir, "roi_diploid.fa")

    tgt_fn = buggy_target if mode == "buggy" else src_calculate_target_pos
    with open(roi_fa, "w") as out:
        for hap, hp in (("hap1", HAP1), ("hap2", HAP2)):
            m = map_paf(seq, hp)
            if not m or m["cigar"] is None:
                print(f"  {hap}: no map"); continue
            tp = tgt_fn(m, voff)
            if tp is None:
                print(f"  {hap}: variant offset not covered"); continue
            # recessive: patch both haps (HOM). dominant: patch hap1 only (HET).
            do_patch = recessive or hap == "hap1"
            s, st = extract_patch(hp, m["contig"], tp, ref, alt, m["strand"], do_patch)
            print(f"  {hap}: contig={m['contig']} strand={m['strand']} target={tp} "
                  f"patched={do_patch} ref-check={st}")
            out.write(f">{hap}_{m['contig']}_{tp}\n{s}\n")

    # simulate paired-end reads
    pref = os.path.join(sample_dir, "reads")
    subprocess.run(["art_illumina", "-i", roi_fa, "-o", pref, "-l", "150",
                    "-f", str(COV), "-p", "-m", "400", "-s", "40", "-q", "-na"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # align to chr15-only hg38 reference, sort, index
    ref_fa = os.path.join(WORK, f"{chrom}.fa")
    if not os.path.exists(ref_fa):
        subprocess.run(["samtools", "faidx", os.path.join(BASE, "hg38.fa"), chrom],
                       check=True, stdout=open(ref_fa, "w"))
        subprocess.run(["samtools", "faidx", ref_fa], check=True)
    bam = os.path.join(sample_dir, "aln.bam")
    p1 = subprocess.Popen(["minimap2", "-ax", "sr", ref_fa,
                           f"{pref}1.fq", f"{pref}2.fq"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-o", bam, "-"], stdin=p1.stdout,
                   check=True, stderr=subprocess.DEVNULL)
    p1.wait()
    subprocess.run(["samtools", "index", bam], check=True)

    # call variants in a window around the truth
    win = f"{chrom}:{pos-60}-{pos+60}"
    vcf = os.path.join(sample_dir, "calls.vcf")
    mp = subprocess.Popen(["bcftools", "mpileup", "-f", ref_fa, "-r", win,
                           "-a", "FORMAT/AD", bam],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(vcf, "w") as vf:
        subprocess.run(["bcftools", "call", "-mv", "-Ov"], stdin=mp.stdout,
                       check=True, stdout=vf, stderr=subprocess.DEVNULL)
    mp.wait()

    # report
    print(f"  --- calls within {win} ---")
    hit = None
    with pysam.VariantFile(vcf) as vf:
        for rec in vf:
            s = rec.samples[0]
            gt = s.get("GT")
            ad = s.get("AD")
            off = rec.pos - pos
            mark = ""
            if rec.pos == pos and ref in (rec.ref,) and alt in rec.alts:
                mark = "  <== EXACT match (pos+allele)"
                hit = ("exact", gt, ad)
            print(f"    {rec.chrom}:{rec.pos} {rec.ref}>{','.join(rec.alts)} "
                  f"GT={gt} AD={ad} offset={off:+d}{mark}")
    # verdict
    print("  --- VERDICT ---")
    if hit:
        gt = hit[1]
        is_hom = gt == (1, 1)
        want = "HOM" if recessive else "HET"
        ok_z = (is_hom == recessive)
        print(f"    position: EXACT ({chrom}:{pos})   "
              f"zygosity: {'HOM' if is_hom else 'HET'} (want {want}) "
              f"{'OK' if ok_z else 'WRONG'}")
        print(f"    => {'PASS' if ok_z else 'FAIL (zygosity)'}")
    else:
        print(f"    position: NO exact call at {chrom}:{pos}  => FAIL (offset)")


if __name__ == "__main__":
    disease = sys.argv[1] if len(sys.argv) > 1 else "tay"
    mode = sys.argv[2] if len(sys.argv) > 2 else "buggy"
    run(disease, mode)
