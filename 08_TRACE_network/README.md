# TRACE network

Construction and scoring of the TRACE regulatory network (V2-2, C-agLINK): a
TF-centred network that links a regulatory SNP to a target gene through three
independent route families.

## Layout

```
08_TRACE_network/
├── call_trace_main.sh      nine copy-paste blocks: six edge tables, Local,
│                           SNP-TF, network, evidence, score, demo
├── TRACE_v2_caller/        all code, with its own README (rules and formats)
└── demo/                   one example file per input type plus a runner
```

## Run

```sh
bash demo/run_demo.sh        # toy dataset, a few seconds
```

For a real run, set the eight input variables in block 01 of `call_trace_main.sh`
and execute the blocks in order. The demo uses the same variable names and the
same output layout, so it doubles as a template.

All six network inputs are built from raw tables by
`TRACE_v2_caller/01_prepare_edge_tables.py`; nothing has to be pre-made. The
network rules, the Route 3 gate and the five final filters are documented in
`TRACE_v2_caller/README.md`.

## Requirements and runtime

Python 3 with duckdb and pandas; no R, no GPU. Large runs need disk space for the
DuckDB checkpoint and spill files.

The bundled demo finishes in seconds. A full run is dominated by the network build
and by converting the evidence parts to text; both scale with the number of
evidence paths.
