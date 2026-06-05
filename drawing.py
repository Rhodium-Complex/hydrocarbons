"""Small graph Laplacian drawing example."""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

def draw_laplacian_example(output_path: str | Path = "output.png") -> None:
    """Draw a simple Laplacian eigenvector plot."""    

    laplacian = np.array(
        [
            [2, -1, 0, 0, -1],
            [-1, 2, -1, 0, 0],
            [0, -1, 2, -1, 0],
            [0, 0, -1, 2, -1],
            [-1, 0, 0, -1, 2],
        ]
    )
    _, eigenvectors = np.linalg.eigh(laplacian)
    eigenvectors = eigenvectors.T

    plt.plot(eigenvectors[1], eigenvectors[2])
    plt.savefig(output_path)


if __name__ == "__main__":
    draw_laplacian_example()
