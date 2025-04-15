import numpy as np
import matplotlib.pyplot as plt

Lap = np.array(
    [[2, -1, 0, 0, -1],
     [-1, 2,-1, 0, 0],
     [ 0,-1, 2,-1, 0],
     [ 0, 0,-1, 2,-1],
     [ -1,-0, 0,-1, 2]])
vals, vecs = np.linalg.eigh(Lap)
vecs = vecs.T
print(vecs)

'''図示する'''

plt.plot(vecs[1], vecs[2])
plt.savefig('output.png')

'''グラフラプラシアンから分子構造を描画する '''
