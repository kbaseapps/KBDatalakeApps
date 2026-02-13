#!/usr/bin/env python3
"""
Enrich pangenome member genomes with coordinates and ontology annotations.

Uses three local query classes:
  - QueryGenomeLocal:     feature coordinates, contig IDs, protein hashes
  - QueryOntologyLocal:   Bakta + KOfam annotations (COG, EC, GO, KO, PFAM, SO, UniRef)
  - (optional) existing RAST .tsv: RAST functional annotations

Writes one TSV per genome: {genome_id}_genome_data.tsv into the same directory
as the FAA files.

Usage (inside the KBase SDK Docker container):
    source /opt/env/berdl_genomes/bin/activate
    python /kb/module/berdl/berdl/bin/enrich_pangenome_features.py \\
        /kb/module/work/shared/chenry/scratch/Acinetobacter_baylyi_ADP1_RAST \\
        RS_GCF_000368685.1
"""
import argparse
import csv
import os
import sys
from pathlib import Path


def _join_set(s):
    """Join a set of annotation values into a semicolon-separated string."""
    if not s:
        return ""
    return ";".join(sorted(s))


def _parse_description(desc):
    """Parse the description string produced by QueryGenomeLocal.get_genome_features().

    Format: '{contig_name} {start} {end} {feature_type} {protein_hash} {cdm_contig_id} {cdm_feature_id} {cdm_protein_id}'
    """
    parts = desc.split()
    if len(parts) >= 5:
        return {
            "contig_id": parts[0],
            "start": parts[1],
            "end": parts[2],
            "feature_type": parts[3],
            "protein_hash": parts[4],
        }
    return {}


def load_rast_tsv(tsv_path):
    """Load an existing 2-column RAST TSV (id, functions) into a dict."""
    rast = {}
    if not os.path.exists(tsv_path):
        return rast
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rast[row["id"]] = row.get("functions", "")
    return rast


def enrich_genome(member_id, query_g, query_o, genome_dir, rast_dict=None):
    """Fetch features + ontology for one genome, return list of row dicts."""
    genome = query_g.get_genome_features(member_id)
    if genome is None:
        print(f"  WARNING: QueryGenomeLocal returned None for {member_id}, skipping")
        return []

    # Bulk ontology lookup
    annotations = query_o.get_protein_ontology_bulk(genome.features)

    rows = []
    for i, feature in enumerate(genome.features):
        meta = _parse_description(feature.description or "")
        annot = annotations.get(i, {})

        rast_function = ""
        if rast_dict and feature.id in rast_dict:
            rast_function = rast_dict[feature.id]

        rows.append({
            "feature_id": feature.id,
            "contig_id": meta.get("contig_id", ""),
            "start": meta.get("start", ""),
            "end": meta.get("end", ""),
            "feature_type": meta.get("feature_type", ""),
            "protein_hash": meta.get("protein_hash", ""),
            "length": len(feature.seq) * 3 if feature.seq else 0,
            "sequence": feature.seq or "",
            "rast_function": rast_function,
            "bakta_function": _join_set(annot.get("bakta_product", set())),
            "cog": _join_set(annot.get("COG", set())),
            "ec": _join_set(annot.get("EC", set())),
            "go": _join_set(annot.get("GO", set())),
            "ko": _join_set(annot.get("KEGG", set())),
            "pfam": _join_set(annot.get("PFAM", set())),
            "so": _join_set(annot.get("SO", set())),
            "uniref_50": _join_set(annot.get("uniref_50", set())),
            "uniref_90": _join_set(annot.get("uniref_90", set())),
            "uniref_100": _join_set(annot.get("uniref_100", set())),
        })

    return rows


TSV_COLUMNS = [
    "feature_id", "contig_id", "start", "end", "feature_type",
    "protein_hash", "length", "sequence",
    "rast_function", "bakta_function",
    "cog", "ec", "go", "ko", "pfam", "so",
    "uniref_50", "uniref_90", "uniref_100",
]


def write_tsv(rows, output_path):
    """Write enriched feature rows to a TSV file."""
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Enrich pangenome member genomes with coordinates and ontology annotations"
    )
    parser.add_argument("output_dir",
                        help="Pipeline output directory (e.g. .../Acinetobacter_baylyi_ADP1_RAST)")
    parser.add_argument("clade_id",
                        help="Pangenome clade ID (e.g. RS_GCF_000368685.1)")
    parser.add_argument("--members", nargs="*", default=None,
                        help="Specific genome IDs to process (default: all members)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    clade_id = args.clade_id
    genome_dir = output_dir / "pangenome" / clade_id / "genome"

    if not genome_dir.exists():
        print(f"ERROR: genome dir does not exist: {genome_dir}")
        sys.exit(1)

    # Discover member genome IDs from existing FAA files
    if args.members:
        member_ids = args.members
    else:
        member_ids = sorted(
            f.stem for f in genome_dir.glob("*.faa")
        )

    if not member_ids:
        print(f"ERROR: no FAA files found in {genome_dir}")
        sys.exit(1)

    print(f"Output dir:  {output_dir}")
    print(f"Clade:       {clade_id}")
    print(f"Genome dir:  {genome_dir}")
    print(f"Members:     {len(member_ids)}")
    print()

    # Initialize query classes
    from berdl.query.query_genome_local import QueryGenomeLocal
    from berdl.query.query_ontology_local import QueryOntologyLocal

    print("Initializing QueryGenomeLocal ...")
    query_g = QueryGenomeLocal()
    print("Initializing QueryOntologyLocal ...")
    query_o = QueryOntologyLocal()
    print()

    for member_id in member_ids:
        print(f"Processing {member_id} ...")

        # Load existing RAST annotations if available
        rast_tsv = genome_dir / f"{member_id}.tsv"
        rast_dict = load_rast_tsv(rast_tsv)
        if rast_dict:
            print(f"  Loaded {len(rast_dict)} RAST annotations from {rast_tsv.name}")

        rows = enrich_genome(member_id, query_g, query_o, genome_dir, rast_dict)

        if not rows:
            print(f"  No features returned, skipping")
            continue

        # Count non-empty annotation fields
        n_bakta = sum(1 for r in rows if r["bakta_function"])
        n_ko = sum(1 for r in rows if r["ko"])
        n_cog = sum(1 for r in rows if r["cog"])
        print(f"  Features: {len(rows)}, bakta: {n_bakta}, ko: {n_ko}, cog: {n_cog}")

        out_path = genome_dir / f"{member_id}_genome_data.tsv"
        write_tsv(rows, out_path)
        print(f"  Wrote: {out_path}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
