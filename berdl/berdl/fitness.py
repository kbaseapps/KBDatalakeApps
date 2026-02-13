import os
import json
import polars as pl
from pathlib import Path
from berdl.hash_seq import ProteinSequence


def map_protein_hash_to_fitness_records(fitness_path='/data/reference_data/phenotype_data/fitness_genomes/'):
    m_to_fitness_feature = {}
    m_to_essentiality = {}  # protein_hash -> set of (genome_id, gene_id)
    essentiality_genome_ids = set()  # genome_ids that have essentiality experiments
    all_fitness_genome_hashes = set()  # ALL protein hashes from fitness genomes
    for genome_id in os.listdir(fitness_path):
        if genome_id.endswith('.json'):
            with open(fitness_path + genome_id, 'r') as fh:
                fitness_genome_data = json.load(fh)
                gid = genome_id[:-5]
                # Check if this genome has essentiality experiments
                genome_has_essentiality = any(
                    'essentiality' in gdata.get('fitness', {})
                    for gdata in fitness_genome_data['genes'].values()
                )
                if genome_has_essentiality:
                    essentiality_genome_ids.add(gid)
                for gene in fitness_genome_data['genes']:
                    _gene_data = fitness_genome_data['genes'][gene]
                    seq = _gene_data.get('protein_sequence')
                    if not seq:
                        continue
                    protein = ProteinSequence(seq)
                    h = protein.hash_value
                    all_fitness_genome_hashes.add(h)
                    _data_fitness = _gene_data.get('fitness')
                    if _data_fitness and len(_data_fitness) > 0:
                        if h not in m_to_fitness_feature:
                            m_to_fitness_feature[h] = {}
                        m_to_fitness_feature[h][(gid, gene)] = {k: (v['fit'], v['t']) for k, v in
                                                                _data_fitness.items() if 'fit' in v}
                        # Track essentiality
                        if 'essentiality' in _data_fitness:
                            if h not in m_to_essentiality:
                                m_to_essentiality[h] = set()
                            m_to_essentiality[h].add((gid, gene))
    return m_to_fitness_feature, m_to_essentiality, essentiality_genome_ids, all_fitness_genome_hashes


def create_genome_fitness_table(input_genome, input_genome_id, m_to_r, r_to_m,
                                m_to_fitness_feature, m_to_essentiality=None,
                                essentiality_genome_ids=None,
                                all_fitness_genome_hashes=None):
    data = {
        'genome_id': [],
        'feature_id': [],
        'feature_h': [],
        'fitness_genome_id': [],
        'fitness_feature_id': [],
        'fitness_feature_h': [],
        'set_id': [],
        'value': [],
    }
    for feature in input_genome.features:
        if feature.seq:
            protein = ProteinSequence(feature.seq)
            h = protein.hash_value
            r = m_to_r[h]
            other_m = r_to_m[r]
            has_any_rows = False
            has_fitness_genome_match = False
            for other_h in other_m:
                if all_fitness_genome_hashes and other_h in all_fitness_genome_hashes:
                    has_fitness_genome_match = True
                genome_gene_fitness = m_to_fitness_feature.get(other_h, dict())
                for (fitness_genome_id, fitness_feature_id), _sets in genome_gene_fitness.items():
                    # Regular fitness records
                    for set_id, (fit, t) in _sets.items():
                        has_any_rows = True
                        data['genome_id'].append(input_genome_id)
                        data['feature_id'].append(feature.id)
                        data['feature_h'].append(h)
                        data['fitness_genome_id'].append(fitness_genome_id)
                        data['fitness_feature_id'].append(fitness_feature_id)
                        data['fitness_feature_h'].append(other_h)
                        data['set_id'].append(set_id)
                        data['value'].append(fit)
                    # Essentiality records for genes in genomes with essentiality data
                    if (m_to_essentiality is not None
                            and essentiality_genome_ids is not None
                            and fitness_genome_id in essentiality_genome_ids):
                        has_any_rows = True
                        is_essential = (
                            other_h in m_to_essentiality
                            and (fitness_genome_id, fitness_feature_id) in m_to_essentiality[other_h]
                        )
                        data['genome_id'].append(input_genome_id)
                        data['feature_id'].append(feature.id)
                        data['feature_h'].append(h)
                        data['fitness_genome_id'].append(fitness_genome_id)
                        data['fitness_feature_id'].append(fitness_feature_id)
                        data['fitness_feature_h'].append(other_h)
                        data['set_id'].append('essentiality')
                        data['value'].append(1.0 if is_essential else 0.0)
            # Marker row: matched to fitness genome cluster but no data
            if has_fitness_genome_match and not has_any_rows:
                data['genome_id'].append(input_genome_id)
                data['feature_id'].append(feature.id)
                data['feature_h'].append(h)
                data['fitness_genome_id'].append('')
                data['fitness_feature_id'].append('')
                data['fitness_feature_h'].append('')
                data['set_id'].append('fitness_genome_match')
                data['value'].append(0.0)

    df = pl.DataFrame(data)
    return df
