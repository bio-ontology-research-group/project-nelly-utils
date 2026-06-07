#!/usr/bin/env python3
"""Confirm the calculate_target_pos coordinate-drift root cause.

For each disease variant: extract the hg38 +/-10kb context, map it to the apr041
hap1/hap2 assemblies with minimap2 -c (CIGAR), then compute the assembly target
position TWO ways:
  (a) buggy linear interpolation  = t_start + (variant_offset - q_start)
  (b) indel-aware CIGAR walk      = true homologous position
delta = (a) - (b) is the coordinate error baked into the simulated genome.
If delta1 != delta2, a homozygous (recessive) patch lands at two different hg38
positions => splits into two heterozygous calls.
"""
import os, sys, json, subprocess, tempfile
import pysam

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HG38 = os.path.join(BASE, "hg38.fa")
HAP1 = os.path.join(BASE, "data", "GCA_050492415.1_apr041.1_v1_genomic.fna")
HAP2 = os.path.join(BASE, "data", "GCA_050492395.1_apr041.2_v1_genomic.fna")
CONFIG = os.path.join(BASE, "config", "diseases.json")
WINDOW = 10000


def get_context(chrom, pos, ref, window=WINDOW):
    with pysam.FastaFile(HG38) as fa:
        if chrom not in fa.references:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt in fa.references:
                chrom = alt
        start = max(0, pos - 1 - window)
        end = pos - 1 + len(ref) + window
        seq = fa.fetch(chrom, start, end)
        return seq, pos - 1 - start  # variant_offset within query


def map_paf(seq, target):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(">query\n" + seq + "\n")
        qp = f.name
    try:
        # -c gives the cg:Z CIGAR tag needed for the indel-aware walk
        out = subprocess.run(
            ["minimap2", "-c", "-x", "asm5", "--secondary=no", target, qp],
            capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(qp)
    for line in out.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        if p[0] == "query":
            cg = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
            return {"q_start": int(p[2]), "q_end": int(p[3]), "strand": p[4],
                    "contig": p[5], "t_start": int(p[7]), "t_end": int(p[8]),
                    "cigar": cg}
    return None


def parse_cigar(cg):
    ops, n = [], ""
    for ch in cg:
        if ch.isdigit():
            n += ch
        else:
            ops.append((int(n), ch))
            n = ""
    return ops


def buggy_target(m, voff):
    """Reproduces scripts/map_and_simulate.py calculate_target_pos (1-based)."""
    if not (m["q_start"] <= voff <= m["q_end"]):
        return None
    off = voff - m["q_start"]
    if m["strand"] == "+":
        return m["t_start"] + off + 1
    return m["t_end"] - off  # = t_end - (off+1) + 1


def cigar_walk_target(m, voff):
    """Indel-aware: walk the CIGAR to the true target pos of query offset voff.
    Returns 1-based target position. CIGAR is in target-forward orientation;
    query walked forward for '+' and backward (from q_end) for '-'."""
    ops = parse_cigar(m["cigar"])
    t = m["t_start"]            # 0-based target cursor
    if m["strand"] == "+":
        q = m["q_start"]
        for length, op in ops:
            if op in ("M", "=", "X"):
                if q + length > voff:           # variant falls inside this block
                    return t + (voff - q) + 1
                q += length; t += length
            elif op == "I":                     # insertion to ref: query advances
                if q + length > voff:            # variant inside an insertion
                    return t + 1                 # collapses onto target cursor
                q += length
            elif op in ("D", "N"):              # deletion: target advances only
                t += length
        return None
    else:
        # reverse: target walks forward from t_start while query walks DOWN from q_end
        q = m["q_end"]
        for length, op in ops:
            if op in ("M", "=", "X"):
                if q - length <= voff:          # variant inside this block
                    return t + (q - voff - 1) + 1
                q -= length; t += length
            elif op == "I":
                if q - length <= voff:
                    return t + 1
                q -= length
            elif op in ("D", "N"):
                t += length
        return None


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    diseases = json.load(open(CONFIG))
    print(f"{'disease':22} {'hap':5} {'strand':6} {'buggy':>11} {'true':>11} {'delta':>6}")
    print("-" * 72)
    for d in diseases:
        if only and only.lower() not in d["name"].lower():
            continue
        v = d["variant"]
        if len(v["ref"]) != 1 or len(v["alt"]) != 1:
            continue  # SNP-only for this check
        try:
            seq, voff = get_context(v["chrom"], v["pos"], v["ref"])
        except Exception as e:
            print(f"{d['name'][:22]:22} context error: {e}")
            continue
        deltas = {}
        for hap, tgt in (("hap1", HAP1), ("hap2", HAP2)):
            m = map_paf(seq, tgt)
            if not m or m["cigar"] is None:
                print(f"{d['name'][:22]:22} {hap:5} NO MAP")
                continue
            b = buggy_target(m, voff)
            tr = cigar_walk_target(m, voff)
            delta = (b - tr) if (b and tr) else None
            deltas[hap] = delta
            print(f"{d['name'][:22]:22} {hap:5} {m['strand']:6} "
                  f"{str(b):>11} {str(tr):>11} {str(delta):>6}")
        if "hap1" in deltas and "hap2" in deltas and deltas["hap1"] is not None and deltas["hap2"] is not None:
            same = deltas["hap1"] == deltas["hap2"]
            print(f"{'':22} -> hap1/hap2 drift {'EQUAL' if same else 'DIFFER'} "
                  f"({deltas['hap1']} vs {deltas['hap2']}) "
                  f"=> {'hom preserved' if same else 'splits into 2 het calls'}")
        print()


if __name__ == "__main__":
    main()
