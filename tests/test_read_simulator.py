import pytest
from src.read_simulator import simulate_reads
import subprocess

def test_simulate_reads_paired_end(mocker):
    """
    Tests the simulate_reads function for paired-end reads, mocking the ART call.
    """
    # Mock shutil.which to always return True
    mocker.patch('shutil.which', return_value=True)
    
    # Mock the subprocess.run call
    mock_run = mocker.patch('subprocess.run', return_value=subprocess.CompletedProcess(args=[], returncode=0, stdout='', stderr=''))

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

    # Check that subprocess.run was called with the correct command
    expected_cmd = [
        'art_illumina',
        '-i', fasta_path,
        '-o', output_prefix,
        '-l', str(read_length),
        '-f', str(coverage),
        '-p',
        '-m', '200',
        '-s', '10'
    ]
    mock_run.assert_called_once_with(expected_cmd, check=True, capture_output=True, text=True)
