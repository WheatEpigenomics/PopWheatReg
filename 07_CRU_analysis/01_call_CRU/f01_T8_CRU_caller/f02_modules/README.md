# Frozen rule modules

The CRU caller is a composition of rule modules that were frozen one by one
during method development. Each file here is the version that produced the
reported calls. No algorithm line was edited; the only change is that the old
`parents[1] / 'fNN_name/0N_script.py'` file loads became relative imports, so the
whole rule set ships as an ordinary Python package.

They must stay together and must keep their names: the modules import each other
at import time, and several expose sub-namespaces (`m11_engine.M`, `.H`) that the
calling chain addresses directly.

| Module | Rule | Imports |
|---|---|---|
| `m01_density_signal.py` | span/endpoint density signal | - |
| `m02_adaptive_background.py` | adaptive background density | - |
| `m03_hard_gap_route.py` | hard-gap B route | - |
| `m04_primary_refine.py` | B1/B2 primary refinement | m03, m02, m01 |
| `m05_density_transitions.py` | agLINK density transitions | m03, m01 |
| `m06_depth_sensitivity.py` | B1-45 depth sensitivity | m05 |
| `m07_changepoint.py` | V6/V7 change point | m06 |
| `m08_baseline_targets.py` | baseline targets and core components | m07, m04 |
| `m09_two_stage_density.py` | V8/V9 two-stage density | m07 |
| `m10_local_subdivision.py` | V11 local subdivision | m09, m02 |
| `m11_engine.py` | V11 engine, exposes `M` and `H` | m10, m08 |
| `m12_block_test.py` | haplotype-block permutation test | - |
| `m13_low_density_recall.py` | T3 low-density recall | m11 (lazy) |
| `m14_support_recall.py` | T4/T5/T6 support recall | m13, m11 (lazy) |
| `m15_island.py` | T7-A island rescue | m14, m01 (lazy) |
| `m16_junctions.py` | shared-junction export | - |
| `m17_block_reference_null.py` | M3 local null generator | - |

`m01` to `m11` build the density profile and cut the genome; `m13` to `m15` add
the low-density recall rounds and the island rescue; `m12`, `m16` and `m17` are
used only by the haplotype-block test.

The historical `main()` entry points and the folder names inside them are
retained for provenance. They are never executed by the CRU caller, which only
uses the rule functions through the namespace handles in `02_call_chain.py`.
