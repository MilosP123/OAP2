"""Run the actual notebook in a fresh Jupyter kernel on synthetic data only."""
import json
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile


def main():
    """Never place synthetic inputs or outputs in the real data directory."""
    project = Path(__file__).resolve().parents[1]
    fixture = runpy.run_path(str(project / 'tests' / 'test_catalog.py'))['make_catalog']
    with tempfile.TemporaryDirectory(prefix='agn-notebook-test-') as temp:
        directory = Path(temp)
        fixture(directory / 'data' / 'group99')
        notebook = directory / 'Untitled7.ipynb'
        shutil.copyfile(project / 'Untitled7.ipynb', notebook)
        result = subprocess.run(
            [sys.executable, '-m', 'pytest', '--nbmake', '--nbmake-timeout=60',
             '-q', str(notebook)], cwd=directory, check=False,
        )
        if result.returncode:
            return result.returncode
        output = directory / 'outputs'
        metadata = json.loads((output / 'catalog_metadata.json').read_text(encoding='utf-8'))
        if metadata['input_kind'] != 'synthetic_test' or metadata['base_sample_count'] != 24:
            raise ValueError('Notebook export does not match the synthetic fixture.')
        print('Fresh-kernel nbmake passed on SYNTHETIC data (not scientific results).')
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
