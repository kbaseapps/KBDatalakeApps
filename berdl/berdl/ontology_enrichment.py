"""
Ontology Term Enrichment from BERDL Data Lake API.

This module provides functions to enrich ontology terms (GO, EC, KEGG, COG, PFAM, SO)
with human-readable labels and definitions from the BERDL ontology source.

Author: Jose P. Faria (jplfaria@gmail.com)
Date: February 2026

=== SUPPORTED ONTOLOGIES ===

- GO (Gene Ontology): GO:XXXXXXX
- EC (Enzyme Commission): EC:X.X.X.X
- KEGG (KO orthologs): KEGG:KXXXXX
- COG (Clusters of Orthologous Groups): COG:COGXXXX, COG:X (categories)
- PFAM (Protein Families): PFAM:PFXXXXX.XX
- SO (Sequence Ontology): SO:XXXXXXX

=== USAGE ===

    from berdl.ontology_enrichment import OntologyEnrichment

    enricher = OntologyEnrichment(token="YOUR_TOKEN")

    # Enrich a list of term IDs
    terms = ["GO:0008150", "EC:1.1.1.1", "KEGG:K00001"]
    enriched = enricher.enrich_terms(terms)
    # Returns: DataFrame with columns [identifier, label, definition]
"""

import time
import requests
import pandas as pd
import re
from typing import List, Dict, Optional, Set


class OntologyEnrichment:
    """
    Enrich ontology terms with labels and definitions from BERDL and other APIs.

    Data sources:
    - BERDL API: GO, EC, SO, PFAM
    - KEGG REST API: KEGG KO (https://rest.kegg.jp/list/ko)
    - Local definitions: COG categories
    """

    BERDL_API_URL = "https://hub.berdl.kbase.us/apis/mcp/delta/tables/query"
    KEGG_API_URL = "https://rest.kegg.jp"

    # API settings
    BATCH_SIZE = 100
    PAGE_LIMIT = 1000
    TIMEOUT = 60

    # Cache for KEGG KO definitions (loaded once)
    _kegg_ko_cache: Dict[str, str] = None

    def __init__(self, token: str):
        """
        Initialize the ontology enrichment client.

        Args:
            token: BERDL API token for KBase authentication.
        """
        if not token:
            raise ValueError("Token is required for BERDL API access.")

        self.token = token
        self.headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
        }

    def _query_berdl(self, query: str) -> List[Dict]:
        """Execute a query against BERDL API with pagination."""
        all_results = []
        offset = 0

        while True:
            payload = {
                "query": query,
                "limit": self.PAGE_LIMIT,
                "offset": offset
            }

            try:
                response = requests.post(
                    self.BERDL_API_URL,
                    headers=self.headers,
                    json=payload,
                    timeout=self.TIMEOUT
                )

                if response.status_code == 200:
                    data = response.json()
                    results = data.get('result', [])
                    all_results.extend(results)

                    if not data.get('pagination', {}).get('has_more', False):
                        break

                    offset += self.PAGE_LIMIT
                else:
                    print(f"BERDL API error {response.status_code}: {response.text[:200]}")
                    break

            except Exception as e:
                print(f"BERDL API exception: {e}")
                break

        return all_results

    def _enrich_from_berdl(self, identifiers: List[str]) -> Dict[str, Dict]:
        """
        Enrich ontology terms using BERDL Data Lake API.
        Works for GO, EC, SO, PFAM.
        """
        if not identifiers:
            return {}

        results = {}

        # Process in batches
        for i in range(0, len(identifiers), self.BATCH_SIZE):
            batch = identifiers[i:i + self.BATCH_SIZE]
            in_clause = ", ".join([f"'{gid}'" for gid in batch])

            query = f"""
                SELECT subject, predicate, value
                FROM kbase_ontology_source.statements
                WHERE subject IN ({in_clause})
                AND predicate IN ('rdfs:label', 'IAO:0000115')
            """

            rows = self._query_berdl(query)

            for row in rows:
                subj = row.get('subject', '')
                pred = row.get('predicate', '')
                val = row.get('value', '')

                if subj not in results:
                    results[subj] = {'label': '', 'definition': ''}

                if pred == 'rdfs:label':
                    results[subj]['label'] = val
                elif pred == 'IAO:0000115':
                    results[subj]['definition'] = val

        return results

    def _load_kegg_ko_list(self) -> Dict[str, str]:
        """
        Load all KEGG KO definitions from KEGG REST API.

        Uses https://rest.kegg.jp/list/ko to get ~26,000 KO definitions in one call.
        Results are cached for subsequent calls.

        Returns:
            Dict mapping KO ID (e.g., 'K00001') to definition string
        """
        if OntologyEnrichment._kegg_ko_cache is not None:
            return OntologyEnrichment._kegg_ko_cache

        print("    Loading KEGG KO definitions from REST API...")
        try:
            url = f"{self.KEGG_API_URL}/list/ko"
            response = requests.get(url, timeout=120)

            if response.status_code == 200:
                ko_dict = {}
                for line in response.text.strip().split('\n'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        # ko:K00001 -> K00001
                        ko_id = parts[0].replace('ko:', '')
                        definition = parts[1]
                        ko_dict[ko_id] = definition

                OntologyEnrichment._kegg_ko_cache = ko_dict
                print(f"    Loaded {len(ko_dict)} KEGG KO definitions")
                return ko_dict
            else:
                print(f"    KEGG API error: {response.status_code}")
                return {}

        except Exception as e:
            print(f"    Error loading KEGG KO list: {e}")
            return {}

    def _enrich_kegg_from_api(self, ko_ids: List[str]) -> Dict[str, Dict]:
        """
        Enrich KEGG KO terms from KEGG REST API.

        Uses the bulk list endpoint for efficiency (~26,000 KOs in one call).
        """
        # Load full KO list (cached after first call)
        ko_definitions = self._load_kegg_ko_list()

        results = {}
        for ko_id in ko_ids:
            # Extract K number (e.g., K00001 from KEGG:K00001)
            k_num = ko_id.replace('KEGG:', '')

            definition = ko_definitions.get(k_num, '')
            if definition:
                # Parse label from definition
                # Format: "E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]"
                # Label is everything before [EC:...] or the whole thing
                label = re.sub(r'\s*\[EC:[^\]]+\]', '', definition).strip()
                results[ko_id] = {'label': label, 'definition': definition}
            else:
                results[ko_id] = {'label': '', 'definition': ''}

        return results

    def _enrich_cog_local(self, cog_ids: List[str], cog_definitions: Dict[str, str] = None) -> Dict[str, Dict]:
        """
        Enrich COG terms from local definitions.

        COG definitions are typically loaded from a local file since they're not
        reliably available via BERDL API.
        """
        results = {}

        # COG category descriptions
        cog_categories = {
            'J': 'Translation, ribosomal structure and biogenesis',
            'A': 'RNA processing and modification',
            'K': 'Transcription',
            'L': 'Replication, recombination and repair',
            'B': 'Chromatin structure and dynamics',
            'D': 'Cell cycle control, cell division, chromosome partitioning',
            'Y': 'Nuclear structure',
            'V': 'Defense mechanisms',
            'T': 'Signal transduction mechanisms',
            'M': 'Cell wall/membrane/envelope biogenesis',
            'N': 'Cell motility',
            'Z': 'Cytoskeleton',
            'W': 'Extracellular structures',
            'U': 'Intracellular trafficking, secretion, and vesicular transport',
            'O': 'Posttranslational modification, protein turnover, chaperones',
            'X': 'Mobilome: prophages, transposons',
            'C': 'Energy production and conversion',
            'G': 'Carbohydrate transport and metabolism',
            'E': 'Amino acid transport and metabolism',
            'F': 'Nucleotide transport and metabolism',
            'H': 'Coenzyme transport and metabolism',
            'I': 'Lipid transport and metabolism',
            'P': 'Inorganic ion transport and metabolism',
            'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
            'R': 'General function prediction only',
            'S': 'Function unknown',
        }

        for cog_id in cog_ids:
            # Handle COG categories (single letter, e.g., COG:J)
            if re.match(r'COG:[A-Z]$', cog_id):
                cat = cog_id.replace('COG:', '')
                results[cog_id] = {
                    'label': cog_categories.get(cat, ''),
                    'definition': ''
                }
            # Handle COG IDs (e.g., COG:COG0001)
            elif cog_definitions and cog_id in cog_definitions:
                results[cog_id] = {
                    'label': cog_definitions[cog_id],
                    'definition': ''
                }
            else:
                results[cog_id] = {'label': '', 'definition': ''}

        return results

    def extract_ontology_terms(self, genome_df: pd.DataFrame) -> Dict[str, Set[str]]:
        """
        Extract ontology term IDs from a genome features DataFrame.

        Looks for columns like 'Annotation:GO', 'Annotation:KO', 'Annotation:EC', etc.

        Returns:
            Dict mapping ontology type to set of term IDs
        """
        terms_by_type = {
            'GO': set(),
            'EC': set(),
            'KEGG': set(),
            'COG': set(),
            'PFAM': set(),
            'SO': set(),
        }

        # Patterns for extracting term IDs
        patterns = {
            'GO': re.compile(r'GO:\d+'),
            'EC': re.compile(r'EC:[\d\.-]+'),
            'KEGG': re.compile(r'(?:KEGG:)?K\d{5}'),
            'COG': re.compile(r'COG:(?:COG\d+|[A-Z])'),
            'PFAM': re.compile(r'(?:PFAM:)?PF\d+(?:\.\d+)?'),
            'SO': re.compile(r'SO:\d+'),
        }

        # Look for annotation columns
        for col in genome_df.columns:
            if not col.startswith('Annotation:'):
                continue

            for _, row in genome_df.iterrows():
                value = str(row.get(col, ''))
                if not value or value == 'nan':
                    continue

                # Try all patterns
                for ont_type, pattern in patterns.items():
                    matches = pattern.findall(value)
                    for match in matches:
                        # Normalize format
                        if ont_type == 'KEGG' and not match.startswith('KEGG:'):
                            match = f'KEGG:{match}'
                        elif ont_type == 'PFAM' and not match.startswith('PFAM:'):
                            match = f'PFAM:{match}'
                        terms_by_type[ont_type].add(match)

        return terms_by_type

    def enrich_terms(self, term_ids: List[str], cog_definitions: Dict[str, str] = None) -> pd.DataFrame:
        """
        Enrich a list of ontology term IDs with labels and definitions.

        Args:
            term_ids: List of ontology term IDs (e.g., ['GO:0008150', 'EC:1.1.1.1'])
            cog_definitions: Optional dict of COG ID -> definition for local COG enrichment

        Returns:
            DataFrame with columns [identifier, label, definition]
        """
        if not term_ids:
            return pd.DataFrame(columns=['identifier', 'label', 'definition'])

        # Group by ontology type
        go_terms = [t for t in term_ids if t.startswith('GO:')]
        ec_terms = [t for t in term_ids if t.startswith('EC:')]
        kegg_terms = [t for t in term_ids if t.startswith('KEGG:')]
        cog_terms = [t for t in term_ids if t.startswith('COG:')]
        pfam_terms = [t for t in term_ids if t.startswith('PFAM:')]
        so_terms = [t for t in term_ids if t.startswith('SO:')]

        all_results = {}

        # Enrich from BERDL (GO, EC, SO, PFAM)
        berdl_terms = go_terms + ec_terms + so_terms + pfam_terms
        if berdl_terms:
            print(f"  Enriching {len(berdl_terms)} terms from BERDL API...")
            berdl_results = self._enrich_from_berdl(berdl_terms)
            all_results.update(berdl_results)

        # Enrich KEGG from KEGG API
        if kegg_terms:
            print(f"  Enriching {len(kegg_terms)} KEGG terms from KEGG API...")
            kegg_results = self._enrich_kegg_from_api(kegg_terms)
            all_results.update(kegg_results)

        # Enrich COG locally
        if cog_terms:
            print(f"  Enriching {len(cog_terms)} COG terms...")
            cog_results = self._enrich_cog_local(cog_terms, cog_definitions)
            all_results.update(cog_results)

        # Build result DataFrame
        rows = []
        for term_id in term_ids:
            info = all_results.get(term_id, {'label': '', 'definition': ''})
            rows.append({
                'identifier': term_id,
                'label': info.get('label', ''),
                'definition': info.get('definition', '')
            })

        return pd.DataFrame(rows)

    def enrich_genome_terms(self, genome_df: pd.DataFrame, cog_definitions: Dict[str, str] = None) -> pd.DataFrame:
        """
        Extract and enrich all ontology terms from a genome features DataFrame.

        Args:
            genome_df: DataFrame with annotation columns (e.g., 'Annotation:GO')
            cog_definitions: Optional dict for COG enrichment

        Returns:
            DataFrame with enriched ontology terms
        """
        # Extract all terms
        terms_by_type = self.extract_ontology_terms(genome_df)

        # Flatten to single list
        all_terms = []
        for terms in terms_by_type.values():
            all_terms.extend(terms)

        if not all_terms:
            print("  No ontology terms found in genome")
            return pd.DataFrame(columns=['identifier', 'label', 'definition'])

        print(f"  Found {len(all_terms)} unique ontology terms")
        return self.enrich_terms(all_terms, cog_definitions)
