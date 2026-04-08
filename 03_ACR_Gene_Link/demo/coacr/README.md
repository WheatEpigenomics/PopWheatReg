## Repository contents
- `04_co_access_analysis.py`  
  Supports downstream co-accessibility-related analysis.

## Demo input

A minimal toy example is provided for the co-accessibility analysis:

- `demo_acr_matrix.txt`: ACR accessibility matrix with ACR IDs as rows and samples as columns.
- `demo_acr.bed`: BED-like file with four columns (`chr`, `start`, `end`, `id`), where `id` matches the row names in the matrix.

## Example command

```bash
python 04_co_access_analysis.py \
  --matrix demo_acr_matrix.txt \
  --bed demo_acr.bed \
  --window 5000 \
  --method pearson \
  --cutoff 0.6 \
  --gpu-mem auto \
  --device cuda:0 \
  --output demo_coaccess_output.tsv
```

## Expected output

The script writes a tab-delimited output file with three columns:

- `id1`
- `id2`
- `correlation`

Only ACR pairs within the specified genomic window and passing the correlation cutoff are retained.
