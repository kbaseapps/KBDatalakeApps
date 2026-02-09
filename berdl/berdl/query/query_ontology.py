from abc import ABC, abstractmethod


class QueryOntologyABC(ABC):

    def get_protein_ontology(self, seq: str):
        """Return all ontologies associated with a protein sequence"""
        raise NotImplementedError
