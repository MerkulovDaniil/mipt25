"""
Experiment 1: VR methods vs SGD on L2-regularized logistic regression.
Demonstrates: GD (fastest) > SVRG ~ SAG (linear) > SGD dec (slow O(1/k)) > SGD const (noise floor).

Key design: moderate condition number (~100-500) so GD converges in ~50 passes,
SVRG/SAG in ~80 passes, SGD dec slowly, SGD const stuck.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from scipy.optimize import minimize

n, d = 2000, 50
mu = 1e-3  # regularization for kappa ~ 50-200

def make_data(seed=42):
    X, y = make_classification(n_samples=n, n_features=d, n_informative=40,
                               n_redundant=10, random_state=seed, class_sep=0.5)
    y = 2 * y - 1
    return X, y

X, y = make_data()

def sigmoid(z):
    return np.where(z >= 0, 1 / (1 + np.exp(-z)), np.exp(z) / (1 + np.exp(z)))

def logistic_loss(w):
    margins = y * (X @ w)
    return np.mean(np.log(1 + np.exp(-np.clip(margins, -500, 500)))) + 0.5 * mu * np.linalg.norm(w)**2

def grad_i(w, i):
    margin = y[i] * (X[i] @ w)
    s = sigmoid(margin) - 1
    return s * y[i] * X[i] + mu * w

def full_grad(w):
    margins = y * (X @ w)
    s = sigmoid(margins) - 1
    return (X.T @ (s * y)) / n + mu * w

# Smoothness constants
eigvals = np.linalg.eigvalsh(X.T @ X / n)
L = eigvals[-1] / 4.0 + mu  # smoothness of full loss
L_max = np.max(np.sum(X**2, axis=1)) / 4.0 + mu  # max per-sample smoothness
kappa = L / mu
print(f"L = {L:.4f}, L_max = {L_max:.4f}, mu = {mu}, kappa = {kappa:.1f}")

# Optimal solution
res = minimize(logistic_loss, np.zeros(d), jac=full_grad, method='L-BFGS-B',
               options={'maxiter': 50000, 'ftol': 1e-15, 'gtol': 1e-14})
f_star = res.fun
print(f"f* = {f_star:.10f}")

n_passes = 100
n_runs = 3

def run_gd(seed=0):
    """Full GD with step 1/L — linear convergence, fastest per-pass."""
    w = np.zeros(d)
    alpha = 1.0 / L
    losses = [logistic_loss(w) - f_star]
    for epoch in range(n_passes):
        g = full_grad(w)
        w = w - alpha * g
        losses.append(logistic_loss(w) - f_star)
    return losses

def run_sgd_const(seed):
    """SGD with constant lr — converges fast then stalls at noise floor ~1e-2."""
    rng = np.random.RandomState(seed)
    w = np.zeros(d)
    # Constant lr: large enough to converge quickly but creates noise floor
    alpha = 1.0 / (10 * L)
    losses = [logistic_loss(w) - f_star]
    for epoch in range(n_passes):
        perm = rng.permutation(n)
        for idx in perm:
            g = grad_i(w, idx)
            w = w - alpha * g
        losses.append(logistic_loss(w) - f_star)
    return losses

def run_sgd_dec(seed):
    """SGD with decreasing lr ~ c/(k+k0) — O(1/k) convergence."""
    rng = np.random.RandomState(seed)
    w = np.zeros(d)
    # Robbins-Monro: alpha_k = c / (k + k0)
    c = 10.0 / L
    k0 = 50.0  # warmup prevents initial instability
    losses = [logistic_loss(w) - f_star]
    for epoch in range(n_passes):
        perm = rng.permutation(n)
        for t_idx, idx in enumerate(perm):
            k = epoch * n + t_idx + 1
            alpha_k = c / (k + k0)
            g = grad_i(w, idx)
            w = w - alpha_k * g
        losses.append(logistic_loss(w) - f_star)
    return losses

def run_sag(seed):
    """SAG — linear convergence, O(np) memory for gradient table."""
    rng = np.random.RandomState(seed)
    w = np.zeros(d)
    alpha = 1.0 / (16 * L_max)
    grad_table = np.zeros((n, d))
    avg_grad = np.zeros(d)
    # Initialize gradient table
    for i in range(n):
        grad_table[i] = grad_i(w, i)
    avg_grad = grad_table.mean(axis=0)

    losses = [logistic_loss(w) - f_star]
    for epoch in range(n_passes):
        perm = rng.permutation(n)
        for idx in perm:
            new_g = grad_i(w, idx)
            avg_grad += (new_g - grad_table[idx]) / n
            grad_table[idx] = new_g
            w = w - alpha * avg_grad
        losses.append(logistic_loss(w) - f_star)
    return losses

def run_svrg(seed):
    """SVRG — linear convergence between SGD and GD.
    Each epoch = 1 full grad + 2n stoch evals = 3 passes.
    """
    rng = np.random.RandomState(seed)
    w = np.zeros(d)
    eta = 1.0 / (3 * L_max)  # step must use per-sample smoothness
    m = 2 * n
    passes_per_epoch = 3
    n_outer = n_passes // passes_per_epoch

    losses = [logistic_loss(w) - f_star]
    for epoch in range(n_outer):
        w_tilde = w.copy()
        g_full = full_grad(w_tilde)
        w_inner = w.copy()
        for t in range(m):
            idx = rng.randint(n)
            v = grad_i(w_inner, idx) - grad_i(w_tilde, idx) + g_full
            w_inner = w_inner - eta * v
        w = w_inner
        val = logistic_loss(w) - f_star
        losses.extend([val] * passes_per_epoch)
    while len(losses) < n_passes + 1:
        losses.append(losses[-1])
    return losses[:n_passes + 1]

# Run
print("Running GD...")
losses_gd = run_gd()

print("Running SGD (const lr)...")
all_sgd = [run_sgd_const(s) for s in range(n_runs)]
losses_sgd = np.median(all_sgd, axis=0)

print("Running SGD (decreasing lr)...")
all_sgd_dec = [run_sgd_dec(s) for s in range(n_runs)]
losses_sgd_dec = np.median(all_sgd_dec, axis=0)

print("Running SAG...")
all_sag = [run_sag(s) for s in range(n_runs)]
losses_sag = np.median(all_sag, axis=0)

print("Running SVRG...")
all_svrg = [run_svrg(s) for s in range(n_runs)]
for lst in all_svrg:
    while len(lst) < n_passes + 1:
        lst.append(lst[-1])
losses_svrg = np.median(all_svrg, axis=0)

# Clip
losses_gd = np.maximum(losses_gd, 1e-16)
losses_sag = np.maximum(losses_sag, 1e-16)
losses_sgd = np.maximum(losses_sgd, 1e-16)
losses_sgd_dec = np.maximum(losses_sgd_dec, 1e-16)
losses_svrg = np.maximum(losses_svrg, 1e-16)

print(f"Final: GD={losses_gd[-1]:.2e}, SAG={losses_sag[-1]:.2e}, SVRG={losses_svrg[-1]:.2e}, "
      f"SGD_dec={losses_sgd_dec[-1]:.2e}, SGD_const={losses_sgd[-1]:.2e}")

# ---- PLOT ----
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'legend.fontsize': 14,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
})

fig, ax = plt.subplots(1, 1, figsize=(10, 6.5))

epochs = np.arange(n_passes + 1)
ax.semilogy(epochs, losses_gd, 'k-', linewidth=2.5, label='GD (полный градиент)', alpha=0.85)
ax.semilogy(epochs, losses_sag, '-', color='#27ae60', linewidth=2.5, label='SAG', alpha=0.85)
ax.semilogy(epochs, losses_svrg, '-', color='#2980b9', linewidth=2.5, label='SVRG', alpha=0.85)
ax.semilogy(epochs, losses_sgd_dec, '--', color='#e67e22', linewidth=2.5,
            label=r'SGD ($\alpha_k \sim 1/k$)', alpha=0.85)
ax.semilogy(epochs, losses_sgd, '-', color='#e74c3c', linewidth=2.5,
            label=r'SGD (пост. шаг)', alpha=0.85)

# Noise floor annotation for SGD const
noise_floor = np.median(losses_sgd[-20:])
if noise_floor > 1e-14:
    ax.axhline(y=noise_floor, color='#e74c3c', linestyle=':', alpha=0.35, linewidth=1)

ax.set_xlabel('Число проходов по данным', fontsize=16)
ax.set_ylabel(r'$f(x^k) - f^*$', fontsize=16)
ax.legend(loc='best', framealpha=0.9, fontsize=14)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, n_passes])
ax.tick_params(labelsize=13)

plt.tight_layout()
outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'files')
plt.savefig(os.path.join(outdir, 'exp_svrg_vs_sgd.pdf'), bbox_inches='tight', dpi=150)
plt.savefig(os.path.join(outdir, 'exp_svrg_vs_sgd.png'), bbox_inches='tight', dpi=150)
print(f"Saved to {outdir}/exp_svrg_vs_sgd.pdf and .png")
