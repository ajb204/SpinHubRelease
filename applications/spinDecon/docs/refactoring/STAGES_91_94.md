# Stages 91-94: domain and workflow consolidation

## Stage 91 - canonical domain identity
`AnalysisMode` and `DatasetTopology` now have one canonical implementation under `domain/`. The historical root modules are compatibility aliases, preserving class identity for external callers. Active project/GUI/processing code uses the canonical topology module.

## Stage 92 - canonical workflow identity
The complete workflow model/status/overview implementations now live only under `workflow/`. Historical root modules are compatibility aliases. The canonical workflow imports the canonical domain `AnalysisMode`, eliminating the duplicate-class identity hazard encountered during the earlier workflow migration.

## Stage 93 - dimension contracts
Dimension guard, labels, peak contract, and viewer contract are canonical under `domain/dimensions/`. Root modules remain compatibility aliases. Active application code now imports the canonical modules.

## Stage 94 - pseudo-axis domain model
`PseudoAxisTable` and related helpers are canonical under `domain/pseudo_axis.py`; the historical root module is a compatibility alias. Active CPMG, decay, and Pseudo2D workspaces use the domain location.

## Regression gate
The full suite is green after each stage. Additional tests explicitly verify that compatibility imports preserve object/class identity, preventing a recurrence of the `AnalysisMode` identity regression.
