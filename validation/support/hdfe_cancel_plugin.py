"""Generate a test-only plugin entry point for deterministic Stata Break tests.

Compile this file in place of src/ctools_plugin.c, using the ordinary plugin
sources and flags. Keep the resulting plugin in a separate test directory.
Production sources never read the injection scalar.
"""
import argparse
from pathlib import Path


def instrument(source):
    original = '    if (SF_poll) SF_poll();\n    return _stata_->stopflag && SW_stopflag;'
    if source.count(original) != 1:
        raise ValueError('Expected exactly one host cancellation callback')
    return source.replace(original, '''    double remaining;
    if (!SF_scal_use("__ctools_test_cancel_after", &remaining) && remaining >= 0) {
        if (remaining <= 1) return 1;
        SF_scal_save("__ctools_test_cancel_after", remaining - 1);
    }
    if (SF_poll) SF_poll();
    return _stata_->stopflag && SW_stopflag;''')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    target = args.output.resolve()
    if target.is_relative_to(root / 'src') or target.is_relative_to(root / 'build'):
        parser.error('Write the test entry point outside src/ and build/')
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(instrument((root / 'src/ctools_plugin.c').read_text()), encoding='utf-8')
