import pytest
from src.read_simulator import simulate_reads
import subprocess

def test_simulate_reads_paired_end(mocker):
    """
    Tests simulate_reads for paired-end mode, mocking the ART call.
    """
    mocker.patch('shutil.which', return_value='/usr/bin/art_illumina')
    mock_run = mocker.patch(
        'subprocess.run',
        return_value=subprocess.CompletedProcess(args=[], returncode=0, stdout='', stderr='')
    )

    fasta_path = "patched.fa"
    output_prefix = "sample1"
    coverage = 30
    read_length = 150

    simulate_reads(
        fasta_path=fasta_path,
        output_prefix=output_prefix,
        coverage=coverage,
        read_length=read_length,
        paired_end=True
    )

    expected_cmd = [
        'art_illumina',
        '-i', fasta_path,
        '-o', output_prefix,
        '-l', str(read_length),
        '-f', str(coverage),
        '-p', '-m', '500', '-s', '50',
    ]
    mock_run.assert_called_once_with(expected_cmd, check=True, capture_output=True, text=True)


def test_simulate_reads_uses_fallback_when_art_missing(tmp_path, mocker):
    """
    Tests that the Python fallback is invoked when art_illumina is not found.
    """
    mocker.patch('shutil.which', return_value=None)
    mock_fallback = mocker.patch('src.read_simulator.simulate_reads_python')

    simulate_reads(
        fasta_path="patched.fa",
        output_prefix=str(tmp_path / "sample"),
        coverage=1,
        read_length=150,
    )

    mock_fallback.assert_called_once()
