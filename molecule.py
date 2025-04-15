import numpy as np

class Molecule:
    """分子構造を表現するクラス。
    
    結合情報（隣接行列）から分子を表現し、固有値を用いた分子識別機能を提供します。
    """
    def __init__(self, bonds: np.ndarray, precision: int = 3):
        """
        分子構造を初期化します。
        
        Args:
            bonds: 結合情報を表す対称な隣接行列
            precision: 固有値の丸め精度（小数点以下の桁数）
        """
        # 入力検証は省略
        
        self.bonds = bonds
        self.eig = np.around(np.sort(np.linalg.eigvalsh(self.bonds)).real, decimals=precision)
        self.fingerprint = tuple(self.eig)
        
    def __hash__(self) -> int:
        """分子構造のハッシュ値を返します。"""
        return hash(self.fingerprint)

    def __eq__(self, other: object) -> bool:
        """2つの分子の等価性を判定します。"""
        if not isinstance(other, Molecule):
            return False
        return self.fingerprint == other.fingerprint
    
    def __len__(self) -> int:
        """分子の原子数を返します。"""
        return len(self.bonds)