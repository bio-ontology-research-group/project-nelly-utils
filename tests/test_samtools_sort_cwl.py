"""Regression tests for the pangenome `samtools-sort.cwl` tool.

Background
----------
On Arvados (2026-06-02) the `samtools-sort` step of `main-vg.cwl` failed with
"Output is missing expected field ... /samtools-sort/aligned_reads_sorted".

Two independent defects produced that failure:

1. **CRAM without a reference.** The command wrote ``samtools sort -O CRAM``
   with no ``--reference``. CRAM is reference-based; on the older ``samtools``
   in the ``coolmaksat/minimap2-samtools`` image this hard-fails
   ("Failed to populate reference for id 0"), so no output file is produced.

2. **Backslash line-continuations inside an interpolated bash block.** CWL
   string interpolation treats ``\\`` as an escape character, so a ``\\`` before
   a newline is *silently stripped* during ``$(...)`` expansion. That collapses
   ``view | sed | sort`` onto one line and bash dies with
   "syntax error near unexpected token '|'".

The fix: pass ``--reference``, run under ``bash -c`` with ``set -euo pipefail``,
write via ``-o``/glob (not stdout), and continue lines with a *trailing pipe*
instead of a backslash.

These tests lock in both fixes. The static guards run anywhere; the end-to-end
test runs only when cwltool + samtools + minimap2 are installed.
"""
import os
import re
import shutil
import subprocess
import random
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Only the deepvariant main-vg.cwl actually wires in a samtools-sort step (it is
# the workflow that failed on Arvados). The vgcall copy of samtools-sort.cwl is
# orphaned — not referenced by its main-vg.cwl or anywhere else — so it is left
# unchanged and intentionally not covered here. Add it back if vgcall ever uses it.
CWL_FILES = [
    os.path.join(
        REPO_ROOT,
        "workflows/pangenome-workflows/GCP/pangenome-GL-wgs-deepvariant-test-pn/samtools-sort.cwl",
    ),
]
# Only check CWL files that actually exist (submodule may not be checked out).
CWL_FILES = [p for p in CWL_FILES if os.path.exists(p)]
requires_cwl = pytest.mark.skipif(
    not CWL_FILES, reason="nelly-workflows submodule not checked out"
)


def _arguments_block(text: str) -> str:
    """Return the `arguments:` block (up to the next top-level key)."""
    m = re.search(r"\narguments:\n(.*?)\n[a-zA-Z]", text, re.DOTALL)
    return m.group(1) if m else ""


@requires_cwl
@pytest.mark.parametrize("cwl", CWL_FILES, ids=lambda p: os.path.basename(os.path.dirname(p)))
def test_cram_output_specifies_reference(cwl):
    """Defect #1: CRAM must be written against an explicit --reference."""
    text = open(cwl).read()
    assert "-O CRAM" in text or "-O, CRAM" in text, "expected CRAM output"
    assert "--reference" in text, (
        "samtools sort -O CRAM without --reference fails on the production image "
        "(no output file) — the reference must be passed explicitly"
    )


@requires_cwl
@pytest.mark.parametrize("cwl", CWL_FILES, ids=lambda p: os.path.basename(os.path.dirname(p)))
def test_no_backslash_continuation_in_pipeline(cwl):
    """Defect #2: no backslash line-continuation inside the interpolated block.

    CWL interpolation strips backslashes, so `\\` + newline breaks the pipeline.
    """
    block = _arguments_block(open(cwl).read())
    assert block, "could not locate arguments block"
    assert "\\\n" not in block, (
        "backslash line-continuation found in the bash pipeline; CWL interpolation "
        "strips it and bash fails with 'syntax error near unexpected token |'. "
        "Use a trailing pipe ('|' at end of line) instead."
    )


@requires_cwl
@pytest.mark.parametrize("cwl", CWL_FILES, ids=lambda p: os.path.basename(os.path.dirname(p)))
def test_uses_pipefail(cwl):
    """A failure in any pipeline stage must fail the step, not be masked."""
    assert "pipefail" in open(cwl).read(), "expected `set -euo pipefail`"


# --------------------------------------------------------------------------- #
# End-to-end run (skipped unless the external tools are available)
# --------------------------------------------------------------------------- #
_TOOLS = ["cwltool", "samtools", "minimap2"]
_missing = [t for t in _TOOLS if shutil.which(t) is None]
requires_tools = pytest.mark.skipif(
    bool(_missing) or not CWL_FILES,
    reason=f"missing tools: {_missing}" if _missing else "submodule not checked out",
)


def _build_fixture(tmp_path):
    """Minimal pangenome-named reference + aligned BAM mimicking the failing input."""
    random.seed(42)
    seq = "".join(random.choice("ACGT") for _ in range(3000))
    wrap = lambda s, w=60: "\n".join(s[i:i + w] for i in range(0, len(s), w))
    ref_pan = tmp_path / "ref_pan.fa"      # contig: GRCh38#0#chr11 (PanSN)
    ref_plain = tmp_path / "ref_plain.fa"  # contig: chr11 (post-sed name)
    reads = tmp_path / "reads.fa"
    ref_pan.write_text(f">GRCh38#0#chr11\n{wrap(seq)}\n")
    ref_plain.write_text(f">chr11\n{wrap(seq)}\n")
    reads.write_text("".join(
        f">r{i}\n{seq[p:p + 150]}\n" for i, p in enumerate(range(100, 2500, 300))
    ))
    subprocess.run(["samtools", "faidx", str(ref_plain)], check=True)
    bam = tmp_path / "aligned.bam"
    p1 = subprocess.Popen(["minimap2", "-a", str(ref_pan), str(reads)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), "-"],
                   stdin=p1.stdout, check=True)
    p1.wait()
    return bam, ref_plain


@requires_tools
@pytest.mark.parametrize("cwl", CWL_FILES, ids=lambda p: os.path.basename(os.path.dirname(p)))
def test_samtools_sort_cwl_end_to_end(cwl, tmp_path):
    """Run the real CWL with cwltool and assert a valid coordinate-sorted CRAM."""
    bam, ref_plain = _build_fixture(tmp_path)

    # Run on host: strip only the DockerRequirement so --no-container is allowed;
    # the command/inputs/outputs are otherwise byte-for-byte the committed tool.
    test_cwl = tmp_path / "samtools-sort.test.cwl"
    test_cwl.write_text("\n".join(
        ln for ln in open(cwl).read().splitlines()
        if "DockerRequirement:" not in ln and "coolmaksat/minimap2-samtools" not in ln
    ))

    job = tmp_path / "job.yml"
    job.write_text(
        f"aligned_reads:\n  class: File\n  path: {bam}\n"
        f"ref:\n  class: File\n  path: {ref_plain}\n"
        f"ref_prefix: GRCh38\n"
    )
    outdir = tmp_path / "out"
    r = subprocess.run(
        ["cwltool", "--no-container", "--outdir", str(outdir), str(test_cwl), str(job)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, f"cwltool failed:\n{r.stderr[-2000:]}"

    cram = outdir / "aligned.sorted.cram"
    assert cram.exists() and cram.stat().st_size > 0, "no CRAM output produced"

    # quickcheck: structurally valid (EOF block present, etc.)
    assert subprocess.run(["samtools", "quickcheck", str(cram)]).returncode == 0

    ref = ["--reference", str(ref_plain)]
    count = subprocess.run(["samtools", "view", *ref, "-c", str(cram)],
                           capture_output=True, text=True).stdout.strip()
    assert count == "8", f"expected 8 records, got {count}"

    header = subprocess.run(["samtools", "view", *ref, "-H", str(cram)],
                            capture_output=True, text=True).stdout
    assert "SO:coordinate" in header, "output is not coordinate-sorted"
    # PanSN prefix must be stripped from contig names.
    assert "SN:chr11" in header and "GRCh38#0#" not in header, "ref_prefix not stripped"
