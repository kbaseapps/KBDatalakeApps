from pathlib import Path
import subprocess


def run_ani(query, library, output_file, threads=20):
    cmd = [
        'skani', 'search', "-t", str(threads),
        "--ql", str(query),
        "-d", str(library),
        "-o", str(output_file),
    ]
    print(' '.join(cmd))
    output = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    return output


def skani_sketch(library: Path, output: Path, t=20):
    if not library.exists() or not library.is_file():
        raise ValueError("invalid library not found/not file")
    cmd = [
        'skani', 'sketch',
        '--fast',
        '-t', str(t),
        '-o', str(output),
        '-l', str(library)
    ]
    output = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    return output
