#!/usr/bin/env python3
"""Merge verified subgroup parquet parts, then remove temporary files after success."""
import argparse
import json
import os
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def _legacy_entries(parts_dir, manifest):
    entries = manifest.get("parts")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Legacy manifest contains no parquet parts.")
    names = [e.get("name") for e in entries]
    expected = [f"part-{i:06d}.parquet" for i in range(1, len(entries) + 1)]
    if names != expected:
        raise ValueError("Legacy part names are missing or out of order.")
    return entries, set(names) | {"COMPLETE.json"}


def _tile_entries(parts_dir, manifest):
    entries = manifest.get("tiles")
    expected_n = manifest.get("expected_tiles")
    if not isinstance(entries, list) or len(entries) != expected_n:
        raise ValueError("Tile manifest is incomplete.")
    run_marker = parts_dir / "RUN.json"
    if not run_marker.is_file():
        raise ValueError(f"Run manifest is missing: {run_marker}")
    run = json.loads(run_marker.read_text())
    n_batches = (run["variant_count"] + run["batch_size"] - 1) // run["batch_size"]
    n_chunks = (run["phenotype_count"] + run["pheno_chunk"] - 1) // run["pheno_chunk"] if run["pheno_chunk"] else 1
    if expected_n != n_batches * n_chunks:
        raise ValueError("Tile count differs from deterministic run layout.")
    allowed = {"RUN.json", "COMPLETE.json"}
    seen = set()
    for index, entry in enumerate(entries):
        batch_id, chunk_id = entry.get("batch_id"), entry.get("chunk_id")
        if not isinstance(batch_id, int) or not isinstance(chunk_id, int):
            raise ValueError("Invalid tile ID in completion manifest.")
        tile = f"tile-g{batch_id:06d}-p{chunk_id:06d}"
        expected_id = (index // n_chunks + 1, index % n_chunks + 1)
        if (batch_id, chunk_id) != expected_id or tile in seen:
            raise ValueError("Tile entries do not match the complete deterministic tile grid.")
        seen.add(tile)
        if entry.get("tile_id") != tile or entry.get("name") != tile + ".parquet" or entry.get("marker") != tile + ".json":
            raise ValueError(f"Invalid tile file names: {tile}")
        marker = parts_dir / entry["marker"]
        if not marker.is_file() or json.loads(marker.read_text()) != entry:
            raise ValueError(f"Tile completion marker is missing or differs: {marker}")
        allowed.update((entry["name"], entry["marker"]))
    return entries, allowed


def merge_parts(out_path, keep_parts=False):
    out_path = Path(out_path)
    parts_dir = Path(f"{out_path}.tmp_parts")
    marker = parts_dir / "COMPLETE.json"
    if not out_path.parent.is_dir():
        raise ValueError(f"Output parent directory does not exist: {out_path.parent}")
    if out_path.exists():
        raise ValueError(f"Output already exists: {out_path}; refusing to overwrite it.")
    if not marker.is_file():
        raise ValueError(f"Incomplete subgroup run: {marker} is missing. "
                         "Temporary parts are preserved; do not merge a partial run.")
    manifest = json.loads(marker.read_text())
    if manifest.get("mode") != "subgroup_omnibus":
        raise ValueError(f"Unexpected subgroup part manifest: {marker}")
    if manifest.get("parts_dir") != str(parts_dir.resolve()):
        raise ValueError(f"Manifest belongs to another parts directory: {marker}")
    version = manifest.get("format_version")
    if version == 1:
        entries, allowed = _legacy_entries(parts_dir, manifest)
    elif version == 2:
        entries, allowed = _tile_entries(parts_dir, manifest)
    else:
        raise ValueError(f"Unsupported subgroup part manifest version: {version}")
    present = {p.name for p in parts_dir.iterdir()}
    if present != allowed:
        raise ValueError(f"Temporary parts directory has missing or unexpected files: {parts_dir}")
    expected_rows = int(manifest["n_rows"])
    if expected_rows < 0 or expected_rows != sum(int(e["rows"]) for e in entries):
        raise ValueError("Manifest row count does not equal the sum of part row counts.")
    schema = None
    for entry in entries:
        part = parts_dir / entry["name"]
        reader = pq.ParquetFile(part)
        if reader.metadata.num_rows != int(entry["rows"]):
            raise ValueError(f"Parquet row count differs from manifest: {part}")
        if schema is None:
            schema = reader.schema_arrow
        elif not reader.schema_arrow.equals(schema):
            reader.close()
            raise ValueError(f"Parquet schema differs between parts: {part}")
        reader.close()
    merging_path = Path(f"{out_path}.merging.tmp")
    if merging_path.exists():
        raise ValueError(f"Incomplete previous merge file exists: {merging_path}")
    writer = None
    try:
        writer = pq.ParquetWriter(merging_path, schema=schema, compression="snappy")
        for entry in entries:
            reader = pq.ParquetFile(parts_dir / entry["name"])
            try:
                for batch in reader.iter_batches(batch_size=65536):
                    writer.write_table(pa.Table.from_batches([batch], schema=schema))
            finally:
                reader.close()
        writer.close()
        writer = None
        actual_rows = pq.ParquetFile(merging_path).metadata.num_rows
        if actual_rows != expected_rows:
            raise ValueError(f"Merged row count {actual_rows:,} != expected {expected_rows:,}; "
                             "temporary parts were preserved.")
        os.replace(merging_path, out_path)
    except Exception:
        if writer is not None:
            writer.close()
        if merging_path.exists():
            merging_path.unlink()
        raise
    if not keep_parts:
        for entry in entries:
            (parts_dir / entry["name"]).unlink()
            if version == 2:
                (parts_dir / entry["marker"]).unlink()
        marker.unlink()
        if version == 2:
            (parts_dir / "RUN.json").unlink()
        for attempt in range(10):
            try:
                parts_dir.rmdir()
                break
            except OSError:
                if attempt == 9 or any(parts_dir.iterdir()):
                    raise
                time.sleep(0.25)  # Lustre may expose the just-deleted entries briefly.
    return expected_rows, len(entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="Final .parquet path used by the subgroup run")
    parser.add_argument("--keep-parts", action="store_true",
                        help="Keep temporary parquet parts after a verified merge")
    args = parser.parse_args()
    try:
        if not args.out.endswith(".parquet"):
            raise ValueError("--out must end in .parquet.")
        n_rows, n_parts = merge_parts(args.out, keep_parts=args.keep_parts)
    except (ValueError, FileNotFoundError, OSError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    print(f"Merged {n_parts:,} parts and {n_rows:,} rows into {args.out}")
    print("Temporary parts kept." if args.keep_parts else "Temporary parts deleted after verification.")


if __name__ == "__main__":
    main()
