import numpy as np


class Molecule:
    """Represent a molecule with its bond matrix and eigenvalue fingerprint."""

    def __init__(self, bonds: np.ndarray, precision: int = 3):
        """Initialize a molecule from a symmetric bond adjacency matrix."""
        self.bonds = bonds
        self.eig = np.around(
            np.sort(np.linalg.eigvalsh(self.bonds)).real,
            decimals=precision,
        )
        self.fingerprint = tuple(self.eig)

    def __hash__(self) -> int:
        """Return the hash of the molecule fingerprint."""
        return hash(self.fingerprint)

    def __eq__(self, other: object) -> bool:
        """Compare molecule fingerprints for equality."""
        if not isinstance(other, Molecule):
            return False
        return self.fingerprint == other.fingerprint

    def __len__(self) -> int:
        """Return the number of atoms in the molecule."""
        return len(self.bonds)
