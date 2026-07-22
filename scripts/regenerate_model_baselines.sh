#!/usr/bin/env bash
# Re-bless the canonical model baselines used by `tests/model_parity.rs`.
#
# Only run this when an intentional algorithmic change is being made — the
# whole point of the baselines is to flag drift that wasn't intended. After
# running, inspect the diff:
#
#   git diff eidolon/test_data/baseline_models/
#
# and confirm the changes are expected before committing.
#
# Usage: ./scripts/regenerate_model_baselines.sh

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
BLESS_BASELINES=1 cargo test --test model_parity
echo
echo "Baselines regenerated. Review the diff before committing:"
echo "  git diff eidolon/test_data/baseline_models/"