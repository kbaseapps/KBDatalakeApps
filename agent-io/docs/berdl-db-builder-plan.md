# BERDL SQLite Database Builder Plan

## Goal
Build a function that compiles an SQLite database from BERDL pipeline output files, matching the schema of the reference E. coli database at `BERDLTablesPrototype/notebooks/data/anl_ecoli/berdl_tables.db`.

## Target Schema (5 tables)

### Table 1: `genome`
| Column | Type | Source |
|--------|------|--------|
| id | VARCHAR(255) PK | `pangenome/{clade}/members.tsv` → genome_id; user genome from `input_params.json` |
| gtdb_taxonomy | VARCHAR(1000) | `members.tsv` → gtdb_taxonomy_id; **MISSING for user genome** |
| ncbi_taxonomy | VARCHAR(1000) | **MISSING** - not in any input file |
| n_contigs | INTEGER | Count contigs from `assembly/*.fna` and `pangenome/{clade}/assembly/*.fna` |
| n_features | INTEGER | Count features from `genome/*.tsv` and `pangenome/{clade}/genome/*.tsv` |

**Missing data:**
- `ncbi_taxonomy` not available for any genome
- `gtdb_taxonomy` not available for the user genome (only for pangenome members)

---

### Table 2: `genome_ani`
| Column | Type | Source |
|--------|------|--------|
| genome1 | VARCHAR(255) PK | `ani/kepangenomes_fast.out` column "Query_name" |
| genome2 | VARCHAR(255) PK | `ani/kepangenomes_fast.out` column "Ref_name" |
| ani | FLOAT | `ani/kepangenomes_fast.out` column "ANI" |
| af1 | FLOAT | `ani/kepangenomes_fast.out` column "Align_fraction_ref" |
| af2 | FLOAT | `ani/kepangenomes_fast.out` column "Align_fraction_query" |
| kind | VARCHAR(255) | Derived from filename: "kepangenomes", "fitness", or "phenotypes" |

**Notes:**
- `fitness_fast.out` and `phenotypes_fast.out` are empty (header only) for ADP1
- Only 1 ANI row will exist (user genome vs clade representative)

---

### Table 3: `genome_features` (most complex)
| Column | Type | Source |
|--------|------|--------|
| id | INTEGER PK | Auto-increment |
| genome_id | VARCHAR(255) | `user_{genome_name}` from `genome/*_kbasedump.tsv` |
| contig_id | VARCHAR(255) | `kbasedump.tsv` → contig column |
| feature_id | VARCHAR(255) | `kbasedump.tsv` → gene_id column |
| length | INTEGER | `abs(end - start)` from kbasedump.tsv |
| start | INTEGER | `kbasedump.tsv` → start |
| end | INTEGER | `kbasedump.tsv` → end |
| strand | VARCHAR(1) | `kbasedump.tsv` → strand |
| sequence | TEXT | `kbasedump.tsv` → protein_translation |
| sequence_hash | VARCHAR(255) | MD5 hash of protein_translation |
| bakta_function | VARCHAR(255) | **MISSING** - no Bakta annotations in ADP1 data |
| rast_function | VARCHAR(255) | `kbasedump.tsv` → functions column |
| cog | VARCHAR(255) | **MISSING** - not in kbasedump.tsv annotations |
| ec | VARCHAR(255) | Parse from SSO annotations in kbasedump.tsv (EC numbers in parentheses) |
| gene_names | VARCHAR(255) | Parse from `kbasedump.tsv` → aliases column (gene field) |
| go | VARCHAR(1000) | **MISSING** - not in kbasedump.tsv annotations |
| ko | VARCHAR(1000) | `kbasedump.tsv` → `Annotation:KO` column |
| pfam | VARCHAR(1000) | `kbasedump.tsv` → `Annotation:PFAM` column |
| so | VARCHAR(1000) | `kbasedump.tsv` → `Annotation:SSO` column |
| uniref_100 | VARCHAR(255) | **MISSING** - no UniRef annotations |
| uniref_90 | VARCHAR(255) | **MISSING** |
| uniref_50 | VARCHAR(255) | **MISSING** |
| pangenome_cluster_id | VARCHAR(255) | Join user features → protein hash → `pangenome_cluster_with_mmseqs.parquet` via mmseqs2 cluster.tsv |
| pangenome_is_core | BOOLEAN | Derived: cluster present in all pangenome genomes |
| psortb | VARCHAR(255) | **MISSING** - no PSORTb data |
| reactions | TEXT | `models/gene_reaction_data.tsv` → reaction column |
| rast_consistency | REAL | **MISSING** - no RAST consistency scores |
| other_rast_annotations | TEXT | **MISSING** |

**Missing data summary for genome_features:**
- `bakta_function` - no Bakta annotations exist
- `cog` - not in KBase genome dump annotations
- `go` - not in KBase genome dump annotations
- `uniref_100/90/50` - no UniRef mapping
- `psortb` - no PSORTb localization data
- `rast_consistency` - no consistency scoring data
- `other_rast_annotations` - no alternate RAST data

**Pangenome cluster assignment challenge:**
The user genome features are NOT in the pangenome parquet. To assign clusters, we need to:
1. Compute protein hashes (SHA-256) for user genome protein sequences
2. Look up those hashes in `master_mmseqs2/master_faa_cluster.tsv` to find representative hashes
3. Match representative hashes back to cluster IDs via the pangenome parquet

---

### Table 4: `missing_functions`
| Column | Type | Source |
|--------|------|--------|
| Reaction | TEXT PK | From `models/genome_reactions.tsv` gapfilling_status != 'none' |
| RAST_function | TEXT | Map reaction → genes → RAST function from genome TSV |
| RichGapfill | INTEGER | 1 if gapfilled under rich media conditions |
| MinimalGapfill | INTEGER | 1 if gapfilled under minimal media conditions |
| PhenotypeGapfill | INTEGER | From `phenotype_tables/genome_phenotypes.tsv` gapfilled_reactions |
| ModuleGapfill | INTEGER | **UNCLEAR** - may come from model data JSON |
| Pangenome | INTEGER | 1 if reaction exists in pangenome member models but not user model |

**Notes:**
- For ADP1, the user model has 0 gapfilled reactions, so this table may be empty or only have pangenome-derived entries
- The `ModuleGapfill` column source is unclear - may need to examine the E. coli pipeline more closely

---

### Table 5: `pan_genome_features`
| Column | Type | Source |
|--------|------|--------|
| id | INTEGER PK | Auto-increment |
| genome_id | VARCHAR(255) | From `pangenome/{clade}/genome/*.tsv` filenames |
| contig_id | VARCHAR(255) | Parse from feature_id (e.g., "NC_005966.1" from "NC_005966.1_1") |
| feature_id | VARCHAR(255) | `pangenome/{clade}/genome/*.tsv` → id column |
| length | INTEGER | Compute from protein sequence length * 3 |
| start | INTEGER | **MISSING** - pangenome genome TSVs only have id,functions |
| end | INTEGER | **MISSING** |
| strand | VARCHAR(1) | **MISSING** |
| sequence | TEXT | From `pangenome/{clade}/genome/*.faa` protein FASTA files |
| sequence_hash | VARCHAR(255) | MD5 of protein sequence |
| cluster_id | VARCHAR(255) | `pangenome_cluster_with_mmseqs.parquet` → cluster_id |
| is_core | BOOLEAN | Derived: cluster present in all pangenome genomes |
| bakta_function | VARCHAR(255) | **MISSING** |
| rast_function | VARCHAR(255) | `pangenome/{clade}/genome/*.tsv` → functions column |
| gene_names | VARCHAR(255) | **MISSING** for pangenome members |
| cog | VARCHAR(1000) | **MISSING** |
| ec | VARCHAR(1000) | Parse from rast_function EC numbers |
| ko | VARCHAR(1000) | **MISSING** |
| pfam | VARCHAR(1000) | **MISSING** |
| go | VARCHAR(1000) | **MISSING** |
| so | VARCHAR(1000) | **MISSING** |
| uniref_100 | VARCHAR(255) | **MISSING** |
| uniref_90 | VARCHAR(255) | **MISSING** |
| uniref_50 | VARCHAR(255) | **MISSING** |

**Missing data summary for pan_genome_features:**
- Only `rast_function` is available from the pangenome genome TSV files
- All other ontology annotations (bakta, cog, ko, pfam, go, so, uniref) are missing for pangenome members
- `start`, `end`, `strand` coordinates are not available for pangenome features
- `length` can be approximated from protein sequence length

---

## Input Directory Structure
```
{output_dir}/
├── input_params.json
├── ani/
│   ├── kepangenomes_fast.out      → genome_ani table
│   ├── fitness_fast.out           → genome_ani table (if non-empty)
│   └── phenotypes_fast.out        → genome_ani table (if non-empty)
├── assembly/
│   └── user_*.fna                 → genome table (contig count)
├── genome/
│   ├── user_*.tsv                 → genome_features (RAST functions)
│   ├── user_*_kbasedump.tsv       → genome_features (main source)
│   └── user_*.faa                 → genome_features (sequences)
├── models/
│   ├── genome_reactions.tsv       → missing_functions, genome_features
│   ├── gene_reaction_data.tsv     → genome_features (reactions column)
│   └── *_data.json                → missing_functions (gapfilling)
├── pangenome/
│   ├── user_to_clade.json         → clade identification
│   └── {clade_id}/
│       ├── members.tsv            → genome table (taxonomy)
│       ├── pangenome_cluster_with_mmseqs.parquet → cluster assignments
│       ├── genome/*.tsv           → pan_genome_features (RAST functions)
│       ├── genome/*.faa           → pan_genome_features (sequences)
│       ├── assembly/*.fna         → genome table (contig count)
│       └── master_mmseqs2/
│           └── master_faa_cluster.tsv → user genome cluster assignment
└── phenotype_tables/
    └── genome_phenotypes.tsv      → missing_functions (phenotype gapfill)
```

## Implementation Strategy

Each table should be implemented as a separate function:
1. `build_genome_table(output_dir, conn)`
2. `build_genome_ani_table(output_dir, conn)`
3. `build_genome_features_table(output_dir, conn)`
4. `build_missing_functions_table(output_dir, conn)`
5. `build_pan_genome_features_table(output_dir, conn)`

Plus a master function: `build_berdl_database(output_dir, db_path)`

## Per-Table Claude Invocation Guide

If running separate Claude instances per table, invoke with:
- "Implement `build_genome_table` following the plan in agent-io/docs/berdl-db-builder-plan.md"
- "Implement `build_genome_ani_table` following the plan in agent-io/docs/berdl-db-builder-plan.md"
- etc.

Each function should be added to `lib/KBDatalakeApps/build_berdl_db.py`.
