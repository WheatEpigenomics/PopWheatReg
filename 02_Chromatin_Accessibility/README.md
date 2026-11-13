# Chromatin accessibility

Peak calling from ATAC-seq and the peak x sample count matrix used by every
downstream step.

## `01_call_peak_and_norm.sh`

Four steps, one block each:

1. **Peak calling.** MACS2 `callpeak` on each library with `--nomodel
   --shift -100 --extsize 200 --scale-to large --keep-dup all -g 16e9`.
2. **Merge and window.** Peaks filtered by IGV are concatenated, sorted, merged,
   kept if at least 200 bp, then cut into 200-bp windows with
   `bedtools makewindows`; windows longer than 150 bp are kept.
3. **Counting.** `bedtools coverage -c` counts fragments per window for every
   library.
4. **Matrix.** `paste` of all `*.filter.count` files gives one row per window and
   one column per library.

Inputs: per-sample spike-normalised fragment BED files (`*_atac_spike.bed`) and
the manually filtered peak BED files (`*filter.peak.bed`).

Outputs: `Merge.filter.bed` (window list) and the merged count matrix.

## Requirements and runtime

MACS2 and bedtools.

No demo is shipped: a run needs MACS2 and the per-library fragment and peak
files. Production scale is about 200 ATAC-seq libraries over 21 chromosomes. The
steps are embarrassingly parallel over libraries.
