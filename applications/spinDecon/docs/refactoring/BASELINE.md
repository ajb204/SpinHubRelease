# Refactor baseline

Baseline established 2026-08-23 before structural refactoring.

## Test invocation
Run from the directory containing the `decon` package:

    PYTHONPATH=. pytest -q decon/tests

## Baseline result
- 297 passed
- 7 failed
- 1 skipped

The seven failures pre-date this refactor and are therefore treated as known baseline failures, not regressions. They concern fit-radius source-path assumptions, pseudo2D diffusion ROI regression expectations, and Slice2D peak-list cleanup expectations.

## Refactor policy
Each structural stage must preserve or improve this baseline. New failures are regressions and block the next stage.
