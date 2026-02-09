from pathlib import Path
import polars as pl
from berdl.query.query_ontology import QueryOntologyABC
from berdl.hash_seq import ProteinSequence
from modelseedpy.core.msgenome import MSFeature


def collect_ontology(_doc):
    ontology = []
    ontology.append(['bakta_product', _doc['product']])
    for db_xref in _doc.get('db_xrefs', []):
        if db_xref.startswith('UniRef:UniRef50'):
            ontology.append(['uniref_50', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('UniRef:UniRef90'):
            ontology.append(['uniref_90', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('UniRef:UniRef100'):
            ontology.append(['uniref_100', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('SO:'):
            ontology.append(['SO', db_xref])
        elif db_xref.startswith('EC:'):
            ontology.append(['EC', db_xref])
        elif db_xref.startswith('KEGG:'):
            ontology.append(['KEGG', db_xref])  # not defined in berld ontology
        elif db_xref.startswith('GO:'):
            ontology.append(['GO', db_xref])
        elif db_xref.startswith('COG:'):
            ontology.append(['COG', db_xref])
        elif db_xref.startswith('PFAM:'):
            ontology.append(['PFAM', db_xref])
        elif db_xref.startswith('UniRef:UniRef100:'):
            ontology.append(['uniref_100', db_xref])
        else:
            ontology.append(['others', db_xref])
    return ontology


class QueryOntologyLocal(QueryOntologyABC):

    def __init__(self, root=Path('/data/reference_data/berdl_db')):
        self.root = root
        self.ldf_annotation_bakta = pl.scan_parquet(root / 'annotation/bakta/*.parquet')
        self.ldf_annotation_kofam = pl.scan_parquet(root / 'annotation/kofam/*.parquet')

    def get_protein_ontology(self, seq: str):
        protein = ProteinSequence(seq)
        h = protein.hash_value

        annotation_bakta = self.ldf_annotation_bakta.filter(pl.col("_id") == h).collect()
        annotation_kofam = self.ldf_annotation_kofam.filter(pl.col("_id") == h).collect()

        return {
            'bakta': annotation_bakta,
            'kofam': annotation_kofam,
        }

    def get_protein_ontology_bulk(self, features_or_seqs: list):
        _all_h = set()
        index_to_h = {}
        for i, f_or_s in enumerate(features_or_seqs):
            seq = f_or_s.seq if isinstance(f_or_s, MSFeature) else f_or_s
            protein = ProteinSequence(seq)
            h = protein.hash_value
            _all_h.add(h)
            index_to_h[i] = h

        collected_annotation = {i: {} for i in range(len(features_or_seqs))}
        # collect bakta ontology
        res = self.ldf_annotation_bakta.filter(pl.col("_id").is_in(_all_h)).collect()
        res = {row['_id']: row for row in res.rows(named=True)}

        for i, h in index_to_h.items():
            feature_annotation_bakta = res.get(h)
            if feature_annotation_bakta:
                feature_ontology = collect_ontology(feature_annotation_bakta)
                for on, v in feature_ontology:
                    if on not in collected_annotation[i]:
                        collected_annotation[i][on] = set()
                    collected_annotation[i][on].add(v)

        # collect kofam ontology
        res = self.ldf_annotation_kofam.filter(pl.col("_id").is_in(_all_h)).collect()
        res = {row['_id']: row['kos'] for row in res.rows(named=True)}

        for i, h in index_to_h.items():
            feature_annotation_kofam = res.get(h, [])
            if feature_annotation_kofam:
                if 'KEGG' not in collected_annotation[i]:
                    collected_annotation[i]['KEGG'] = set()
                for ko in feature_annotation_kofam:
                    collected_annotation[i]['KEGG'].add(ko)

        return collected_annotation
