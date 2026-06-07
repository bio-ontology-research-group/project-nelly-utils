#!/usr/bin/env python3
"""Scale the end-to-end validation across every SNP disease in the panel.

Runs validation/roundtrip.py in the given mode for each disease and prints a
single pass/fail table. Each disease: simulate (haplotype round-trip) -> align
to hg38 -> call -> compare recovered variant to the planted truth.
"""
import os, sys, json, subprocess, re

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RT = os.path.join(BASE, "validation", "roundtrip.py")
CONFIG = os.path.join(BASE, "config", "diseases.json")

mode = sys.argv[1] if len(sys.argv) > 1 else "fixed"
diseases = json.load(open(CONFIG))

rows = []
for d in diseases:
    v = d["variant"]
    if len(v["ref"]) != 1 or len(v["alt"]) != 1:
        rows.append((d["name"], "SKIP (not a SNP — repeat/indel)", ""))
        continue
    key = d["name"].split()[0].lower()
    print(f"# running {d['name']} ({mode}) ...", flush=True)
    out = subprocess.run([sys.executable, RT, key, mode],
                         capture_output=True, text=True).stdout
    verdict = "NO VERDICT"
    detail = ""
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("=> "):
            verdict = line[3:]
        if line.startswith("position:"):
            detail = line
    rows.append((d["name"], verdict, detail))

print("\n" + "=" * 78)
print(f"END-TO-END SIMULATOR VALIDATION  (mode={mode})")
print("=" * 78)
npass = nfail = 0
for name, verdict, detail in rows:
    print(f"{name:32} {verdict}")
    if detail:
        print(f"{'':32}   {detail}")
    if verdict.startswith("PASS"):
        npass += 1
    elif verdict.startswith("FAIL"):
        nfail += 1
print("-" * 78)
print(f"PASS={npass}  FAIL={nfail}  (SNP diseases tested)")
