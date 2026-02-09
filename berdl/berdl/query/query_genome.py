from abc import ABC, abstractmethod
from typing import Iterable, Protocol, Optional, Any, Union


class QueryGenomeABC(ABC):

    @abstractmethod
    def get_genome_features(self, genome_id: str) -> Any:
        """Return all genome features."""
        raise NotImplementedError

    @abstractmethod
    def get_genome_contigs(self, genome_id: str) -> Any:
        """Return all contigs."""
        raise NotImplementedError
