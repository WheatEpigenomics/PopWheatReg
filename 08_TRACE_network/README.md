# TRACE (Tracing Regulatory Attribution across Chromatin and Expression) network

Construction and scoring of the TRACE regulatory network (V2-2, C-agLINK), a TF-centred framework that connects regulatory SNPs to target genes through three complementary regulatory routes, integrating genetic associations, chromatin accessibility, ACR–gene links and TF binding information.

## TRACE framework

<p align="center">
  <img width="750" alt="TRACE framework" src="https://github.com/user-attachments/assets/364a0475-97db-4cdb-9169-c27df5f47811" />
</p>

## Workflow overview

The TRACE pipeline consists of four main stages:

1. **Edge preparation:** Standardize eQTL, caQTL, agLINK, coACR and TF-binding data into network edge tables, together with SNP functional annotations.
2. **SNP–TF association:** Identify candidate SNP–TF relationships using local regulatory units (L1/L2) and proximal SnpEff annotations.
3. **Regulatory network construction:** Integrate three regulatory routes to establish SNP–TF–gene relationships:
   - **Trans-eQTL route:** SNP–TF associations supported by trans-eQTL evidence.
   - **Trans-caQTL route:** SNP–ACR–gene regulatory relationships supported by trans-caQTL evidence.
   - **TF-binding route:** SNP–TF associations connected to target genes through TF binding and ACR–gene regulatory links.
4. **Evidence filtering and scoring:** Apply route-specific criteria and five final filters, aggregate supporting evidence, and score TF–gene relationships by the number of distinct regulatory routes.

Detailed network rules, including the C-agLINK SNP–TF gate, distance criteria and evidence filtering, are provided in [TRACE_v2_caller/README.md](TRACE_v2_caller/README.md).

## Layout

```text
08_TRACE_network/
├── call_trace_main.sh      Main pipeline with eight numbered sections
├── TRACE_v2_caller/        Network construction and scoring scripts
│   └── README.md           Detailed rules, input formats and outputs
└── demo/                   Example inputs and demonstration runner
```

## Run

Run the bundled demonstration from the `08_TRACE_network/` directory:

```sh
bash demo/run_demo.sh
```

For a full analysis, configure the eight input variables in section 01 of `call_trace_main.sh` and execute the subsequent sections in order.

The pipeline prepares standardized network inputs, constructs local regulatory units and SNP–TF relationships, assembles the regulatory network, filters and consolidates evidence, and calculates network scores.

The demo follows the same workflow and provides example files for each input type.

## Outputs

The primary outputs include:

- **TF–gene network:** Regulatory relationships with supporting routes and evidence counts (`03_TF_network.network.tsv`).
- **Regulatory evidence:** SNP-level supporting evidence for each inferred relationship (`04_TF_network.evidence.tsv.gz`).
- **Network scores:** TF–gene relationships ranked by the number of distinct supporting regulatory routes (`06_TF_network_main_scored_simple.tsv`).

Additional quality-control tables report SNP–TF entrance decisions and filtering results.

## Requirements and runtime

Python 3.8+ with `duckdb` and `pandas`; no R or GPU required. Large-scale analyses require sufficient disk space for DuckDB checkpoints and temporary files.

The bundled demo completes in seconds. Full-scale runtime is primarily determined by network construction and evidence conversion, depending on the number of regulatory associations and supporting paths.
