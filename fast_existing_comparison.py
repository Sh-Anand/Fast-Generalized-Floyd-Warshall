import numpy as np
from hwcounter import count, count_end
import sys
import matplotlib.pyplot as plt

#source : https://github.com/mattianeroni/Fast-Floyd-Warshall
def floyd_warshall(A):
    """
    A is the matrix of distances between nodes.
    """
    n,m = A.shape
    I = np.identity(n)
    A[I == 1] = 0 # diagonal elements should be zero
    for i in range(n):
        A = np.minimum(A, A[i,:] + A[:,i])

def measure_fast(n):
    cycles = []
    A = np.random.random((n, n))
    for i in range(5):
        start = count()
        floyd_warshall(A)
        end = count_end() - start
        cycles.append(end)
    return np.mean(cycles)

x = [8, 24, 72, 216, 648, 1944]
ours = [1.475079, 3.189038, 4.230244, 4.094568, 3.704055, 3.592653]

fast = []

for n in x:
    cycles = measure_fast(n)
    fast.append(2*(n**3)/cycles)

plt.plot(x, ours, label="tiled vectorized", marker='o')
# plt.plot(x, y[fwi][1], label="basic optimizations", marker='o')
plt.plot(x, fast, label="fast python implementation", marker='o')

scalar_peak_performance = 2.48

# plt.axhline(y=vect_peak_performance, label="vectorized peak performance", color='r', linestyle='--')
plt.axhline(y=scalar_peak_performance, label="scalar peak performance", color='y', linestyle='--')

# Gray background lines
for i in range(0, 11):
    plt.axhspan(i, i+.5, facecolor='0.1', alpha=0.1)

plt.suptitle("min plus comparison with fast existing implementation", fontweight="bold", fontsize=12, x=0.31)
plt.legend()
plt.xscale('log', base=3)
plt.ylim([0.0, 5.5])
plt.xlabel("N")
plt.ylabel("flops/cycle")

plt.savefig("./plots/" + "fast_existing_comparison_performance_plot.pdf", bbox_inches="tight")
plt.cla()
