#!/usr/bin/env python3
"""Small stdlib checks of qualification report failure semantics."""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name('qualify_summary.py')


class SummaryTests(unittest.TestCase):
    def run_summary(self, data, strict=False):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'fixtures.tsv'
            log.write_text(data)
            return subprocess.run(
                [sys.executable, '-I', str(SCRIPT), str(log)]
                + (['--require-all'] if strict else []),
                text=True, capture_output=True, check=False)

    def test_counts_and_strict_skips(self):
        data = 'a\tone\texecuted\tinput\nb\tone\tskipped\toracle\nc\ttwo\texecuted\tbinary\n'
        result = self.run_summary(data)
        self.assertEqual(result.returncode, 0)
        self.assertIn('one: executed=1 skipped=1', result.stdout)
        self.assertIn('TOTAL: executed=2 skipped=1', result.stdout)
        self.assertIn('SKIPPED one::b: oracle', result.stdout)
        self.assertEqual(self.run_summary(data, True).returncode, 1)

    def test_empty_and_invalid(self):
        self.assertEqual(self.run_summary('').returncode, 1)
        self.assertEqual(self.run_summary('a\tone\tskipped\tx\n').returncode, 1)
        self.assertEqual(self.run_summary('torn record').returncode, 2)
        self.assertEqual(self.run_summary('a\tone\texecuted\tx').returncode, 2)
        self.assertEqual(self.run_summary('a\tone\tpassed\tx\n').returncode, 2)

    def test_all_executed(self):
        self.assertEqual(self.run_summary('a\tone\texecuted\tx\n', True).returncode, 0)


    def test_entrypoint_propagates_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake_bin = root / 'bin'
            fake_bin.mkdir()
            cargo = fake_bin / 'cargo'
            cargo.write_text("""#!/usr/bin/env bash
case "$1" in
  --version) echo 'cargo mock';;
  build) exit "${MOCK_BUILD_EXIT:-0}";;
  test) printf 'example\\tmock\\t%s\\tfixture\\n' "${MOCK_DECISION:-executed}" >> "$NOMICS_FIXTURE_LOG"; exit "${MOCK_TEST_EXIT:-0}";;
esac
""")
            cargo.chmod(0o755)
            cases = [('pass', 0, 0, 'executed', 0),
                     ('build-fail', 1, 0, 'executed', 1),
                     ('test-fail', 0, 1, 'executed', 1),
                     ('skip', 0, 0, 'skipped', 1)]
            for name, build, tests, decision, expected in cases:
                env = os.environ.copy()
                env.update(PATH=str(fake_bin) + os.pathsep + env['PATH'],
                           MOCK_BUILD_EXIT=str(build), MOCK_TEST_EXIT=str(tests),
                           MOCK_DECISION=decision)
                out = root / name
                result = subprocess.run(
                    ['bash', str(SCRIPT.with_name('qualify.sh')), str(out)],
                    env=env, text=True, capture_output=True, check=False)
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertTrue((out / 'report.txt').is_file())


if __name__ == '__main__':
    unittest.main()
