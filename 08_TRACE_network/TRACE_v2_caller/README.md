
# TRACE (Tracing Regulatory Attribution across Chromatin and Expression) network

Construction of the TRACE regulatory network: a TF-centred network that links a
regulatory SNP to a target gene through three independent route families.
`call_trace_main.sh` one directory up holds the runnable blocks; this directory
holds the code.

The V2-2 rules are the final version used for the csv2 network. Only code is
shipped; input files and results are not.

## TRACE framework
<p align="center">
  <img width="750" alt="Trace_frame" src="https://github.com/user-attachments/assets/364a0475-97db-4cdb-9169-c27df5f47811" />
</p>

## Requirements

- Python 3.8+ with `duckdb` and `pandas`
- Inputs prepared by the user in the formats below

Large runs write DuckDB checkpoints and spill files, so reserve disk space
accordingly and reuse the same `--checkpoint-db` with `--resume` after an
interruption.

## Input formats

| Raw input | Format | Converted to |
|---|---|---|
| eQTL | 9-column TSV, no header: `chr  start  end  SNP  gene  P  cis/tra  beta  beta_se` | SNP-Gene edge |
| caQTL | TSV with header, must contain `variant_id, phenotype_id, pvalue, cis_trans, beta, beta_se, methods, p_source, locus_id` | SNP-ACR edge |
| agLINK | 5-column TSV, no header: `gene  ACR  value  p_or_NA  aglink|calink|dist` | ACR-Gene edge |
| coACR | TSV: `id1  id2  correlation`, header line optional | ACR-ACR edge |
| TF binding | TSV, header line optional: `TF  motif  ACR  PCC  effect` | TF-ACR edge |
| SnpEff ANN | TSV with `ID` plus `ANN.GENE / ANN.EFFECT / ANN.IMPACT` (or the `ANN[0].*` variant) | ANN + `canonical_gene_id` |
| gene BED | BED6: `chr  start  end  gene_id  score  strand` | used directly |
| TF list | TSV with a `TF` column | used directly |

For agLINK, `aglink`/`calink` rows become `agLINK` edges with `score = value` and
`meta = PCC=value`; `dist` rows become `distance` edges with `score = NA` and
`meta = dist=value`. For binding, the score is `abs(PCC)` while `meta` keeps the
signed value: `PCC=..;effect=..;motif=..`.

Conventions:

- SNP IDs are `chr:pos`, 1-based (`chr_pos` and `chr-pos` are also accepted).
- ACR IDs are `name:chr_start_end` with 0-based half-open BED coordinates.
- `cis`/`tra` labels are trusted as given, but a same-chromosome `tra` row must
  be more than 1 Mb away; the converter raises an error otherwise instead of
  silently relabelling it.
- caQTL must be the **source-annotated** table: the `beta`/`beta_se` pair is the
  representative effect matching `p_source`. A PW record has no single signed
  effect and keeps `NA`.
- The edge converters refuse to overwrite an existing output; delete the file or
  pick a new name.
- Binding keeps only `abs(PCC) >= --min-abs-pcc` (0.35 by default) and drops ACRs
  on `Un`, matching the published table; `--keep-un` disables the latter.
- The full-precision `abs(PCC)` is written to `score`. The historical table had
  six significant digits there because of an awk number-formatting artefact; the
  signed `PCC` inside `meta`, which is what the network reads, is identical.

## Steps

```
01_prepare_edge_tables.py   edges            eQTL/caQTL table -> 7-column SNP-Gene / SNP-ACR edge
                            aglink           agLINK+distance table -> ACR-Gene edge
                            coacr            co-accessibility table -> ACR-ACR edge
                            binding          TF-motif-ACR PCC table -> TF-ACR edge
                            normalize-ann    append canonical_gene_id to the ANN table
02_build_local_units.py                      L1/L2 local units from the four edge tables
03_build_snp_tf_edges.py                     SNP-TF edges from Local units + proximal ANN
05_build_caglink_network.py                  V2-2 network build (rules below)
06_finalize_evidence.py                      evidence parts -> one gzipped TSV
07_score_network.py                          simple Route-count score
```

`04_network_rules.py` is a library, not an entry point: it holds the route SQL,
the five final filters and the network assembly used by the builder.

## Rules

### Three routes

| Route | Path |
|---|---|
| `eqtl_trans` | SNP —(trans-eQTL)→ TF gene; TF → target gene |
| `caqtl_trans` | SNP —(trans-caQTL)→ ACR; ACR → target gene |
| `binding` | SNP–TF pair (gate below); TF —(binding)→ ACR; ACR —(direct agLINK)→ target gene |

The ACR hub prefers a direct agLINK over a distance edge for the same ACR–gene
pair, and only falls back to a coACR hop when no direct pair exists.

### Route 3 gate (the C-agLINK rule)

A SNP–TF pair is kept when

```
has_local_unit
  OR ( has_snpeff_prox
       AND within_5kb
       AND ( has_cis_eqtl OR n_aglink_cis_acrs > 0 ) )
```

- `within_5kb`: distance from the SNP to the nearest base of the TF gene body is
  at most 5 kb on the same chromosome.
- `has_cis_eqtl`: a cis-eQTL record exists for the same SNP and the same TF gene.
- `n_aglink_cis_acrs`: at least one ACR is reachable from that SNP by cis-caQTL
  and is linked to the TF gene by a direct agLINK, with
  `max(SNP, ACR midpoint, TSS) - min(...) <= 500 kb`.

All three conditions are required together; a proximal annotation alone is not
enough, and a cis-caQTL chain alone is not enough either.

### Final filters

Each SNP-level evidence row must pass all five tests before aggregation:

| Flag | Rule |
|---|---|
| `pass_coordinates` | TF, gene, SNP coordinates resolvable; ACR resolvable when present |
| `pass_no_self` | TF is not the target gene |
| `pass_tf_gene` | same-chromosome TF–gene **midpoint** distance > 500 kb |
| `pass_tf_acr` | same-chromosome TF–ACR **midpoint** distance > 500 kb (passes when there is no ACR) |
| `pass_snp_gene` | same-chromosome SNP to nearest gene-body base > 500 kb; a different chromosome passes |

Rows are filtered at SNP level first, then aggregated to SNP–TF–gene and finally
to TF–gene, so a TF–gene relation keeps the support of its passing SNPs only.

Three distance conventions coexist and must not be mixed: the upstream 1 Mb
`cis`/`tra` label, the Local span (SNP / ACR floor midpoint / strand-aware TSS,
<= 500 kb) and the final filters above (midpoints for TF pairs, gene-body
interval for SNP–gene).

### Score

`score` is the number of distinct route labels behind one TF–gene relation, not
the number of evidence paths. The table is sorted by `score` descending and then
`n_support` descending; `net_id` is generated once by the builder and never
renumbered.

## Outputs

| Path | Content |
|---|---|
| `03_TF_network.network.tsv` | one row per TF–gene relation, native `TFN000000001` ids |
| `03_TF_network.evidence.part*.parquet` | full evidence paths |
| `03_TF_network.filter_qc.tsv` | per-rule and per-route losses |
| `01_Caglink_SNP_TF_entrances.tsv`, `02_Caglink_entrance_QC.tsv` | Route 3 gate decisions and counts |
| `04_TF_network.evidence.tsv.gz` | evidence after conversion, gzip CRC and SHA256 verified |
| `06_TF_network_main_scored_simple.tsv` | scored main table |

## Demo

```sh
bash demo/run_demo.sh
```

`demo/data/` holds one tiny example of every input file, so the formats above
can be read literally. The demo runner uses the same variable names, the same
seven blocks and the same output layout as `call_trace_main.sh`; only the paths
differ, so it doubles as a template for a real run.

| File | Shows |
|---|---|
| `genes.bed` | two TF genes and two target genes, both strands |
| `TF_subfamilies.txt` | the TF column |
| `snpeff_allANN.tsv` | the ANN columns, one proximal record per SNP |
| `eqtl_fil_beta.bed` | cis and trans rows for both SNPs |
| `caqtl_source_beta.tsv.gz` | cis row plus trans rows, one of them a PW record with `methods=Global;sg1`; gzip input |
| `aglink_merge_distance.txt` | `aglink` and `dist` rows for the same gene |
| `coacces_acr.tsv` | one co-accessibility pair |
| `tf_acr_pcc.tsv` | two signals above the cutoff, one below it |

The demo writes 6 evidence rows, 3 SNP-TF-gene relations and 3 TF-gene
relations covering all three routes, and finishes in a few seconds:

```
net_id        snp        tf                   gene                 routes                         score
TFN000000003  1A:200500  TraesCS1A03G0000200  TraesCS1A03G0000400  binding|caqtl_trans|eqtl_trans  3
TFN000000001  1A:102000  TraesCS1A03G0000100  TraesCS1A03G0000300  binding|eqtl_trans              2
TFN000000002  1A:102000  TraesCS1A03G0000100  TraesCS1A03G0000400  caqtl_trans                     1
```

The first row also shows the two SNP-TF sources merged (`local_unit|tf_prox_snpeff`)
and the pw provenance surviving into `qtl_methods` / `qtl_p_source`.
