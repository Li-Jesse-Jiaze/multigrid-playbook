import numpy as np
import matplotlib.pyplot as plt

# 1D multigrid error smoothing process using weighted Jacobi relaxation
# --------------------------------------------------------------
# Grid size and spacing
N = 80
h = 1.0 / (N + 1)
# 1D Poisson operator (Dirichlet BCs)
A = (np.diag(2*np.ones(N)) +
     np.diag(-1*np.ones(N-1), k=1) +
     np.diag(-1*np.ones(N-1), k=-1)) / h**2

# Weighted Jacobi parameter (optimal ~2/3 for Poisson)
omega = 2.0 / 3.0
# Inverse of diagonal of A
D_inv = np.diag(1.0 / np.diag(A))
# Smoothing operator S = I - omega * D^{-1} A
S = np.eye(N) - omega * (D_inv @ A)

# Initial error: random high-frequency content
e0 = np.random.randn(N)

n_sweeps = 10  # number of smoothing sweeps

# Collect error vectors for each sweep
errors = [e0.copy()]
e = e0.copy()
for sweep in range(n_sweeps):
    e = S @ e
    errors.append(e.copy())

# Spatial grid for plotting
x = np.linspace(h, 1-h, N)

plt.figure(figsize=(12, 6))
# Plot error at each sweep
for i, err in enumerate(errors):
    alpha = 0.7 if 0 < i < n_sweeps else 1.0
    lw = 1.5
    plt.plot(x, err, alpha=alpha, linewidth=lw)

# Remove axes, grid, and legend
plt.axis('off')

plt.tight_layout()
plt.savefig("../slides/fig/cover.svg", transparent=True)
# plt.show()