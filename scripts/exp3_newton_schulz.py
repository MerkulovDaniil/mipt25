"""
Experiment 3: Newton-Schulz approximation quality + effect on optimization
Left: NS convergence to polar factor for different condition numbers
Right: How many NS iterations are needed for Muon? (with EMA smoothing)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(42)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# --- Left panel: NS convergence ---
ax = axes[0]
m, n_dim = 30, 20

cond_numbers = [1, 10, 100, 1000]
colors_left = ['#2ecc71', '#3498db', '#e67e22', '#e74c3c']
markers = ['o', 's', '^', 'D']

for kappa, color, marker in zip(cond_numbers, colors_left, markers):
    rng = np.random.RandomState(42)
    U_ = np.linalg.qr(rng.randn(m, m))[0][:, :n_dim]
    V_ = np.linalg.qr(rng.randn(n_dim, n_dim))[0]
    s = np.linspace(1, kappa, n_dim)
    G = U_ @ np.diag(s) @ V_.T

    # True polar factor
    U_true, _, Vt_true = np.linalg.svd(G, full_matrices=False)
    P_true = U_true @ Vt_true

    # Newton-Schulz
    a, b, c = (3.4445, -4.7750, 2.0315)
    X = G / np.linalg.norm(G, 'fro')
    X = X / (max(G.shape)**0.5)

    errors = []
    for it in range(12):
        scale = np.linalg.norm(P_true, 'fro') / (np.linalg.norm(X, 'fro') + 1e-12)
        err = np.linalg.norm(X * scale - P_true, 'fro') / np.linalg.norm(P_true, 'fro')
        errors.append(max(err, 1e-16))
        A_ = X @ X.T
        B_ = b * A_ + c * A_ @ A_
        X = a * X + B_ @ X

    ax.semilogy(range(12), errors, f'-{marker}', color=color, linewidth=2.5, markersize=7,
                label=f'$\\kappa$ = {kappa}', alpha=0.85)

ax.set_xlabel('Newton-Schulz iteration', fontsize=14)
ax.set_ylabel(r'Relative error $\|X_k - UV^\top\| / \|UV^\top\|$', fontsize=13)
ax.set_title('Convergence to polar factor', fontsize=14)
ax.legend(fontsize=12, loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 11])
ax.tick_params(labelsize=12)
ax.axhline(y=0.01, color='gray', linestyle=':', alpha=0.4)
ax.text(8, 0.016, '1% error', color='gray', fontsize=10, alpha=0.6)

# Annotation: 5 iterations is enough
ax.axvline(x=5, color='#9b59b6', linestyle='--', alpha=0.5, linewidth=1.5)
ax.text(5.2, 2e-8, 'Muon default\n(5 iters)', color='#9b59b6', fontsize=10, alpha=0.7)

# --- Right panel: effect on optimization (with EMA smoothing) ---
ax = axes[1]

m_opt, d_opt, p_opt = 200, 20, 10
rng = np.random.RandomState(123)
U_A = np.linalg.qr(rng.randn(m_opt, m_opt))[0][:, :d_opt]
V_A = np.linalg.qr(rng.randn(d_opt, d_opt))[0]
s_A = np.logspace(0, 2, d_opt)
A_opt = U_A @ np.diag(s_A) @ V_A.T
B_opt = rng.randn(m_opt, p_opt) * 0.1
n_samples = m_opt

def loss(W):
    return 0.5 * np.linalg.norm(A_opt @ W - B_opt, 'fro')**2 / n_samples

W_star = np.linalg.lstsq(A_opt, B_opt, rcond=None)[0]
f_star = loss(W_star)

def newton_schulz(G, num_iters):
    a, b, c = (3.4445, -4.7750, 2.0315)
    X = G / (np.linalg.norm(G, 'fro') + 1e-12)
    X = X / (max(G.shape)**0.5)
    for _ in range(num_iters):
        A_ = X @ X.T
        B_ = b * A_ + c * A_ @ A_
        X = a * X + B_ @ X
    return X

def stoch_grad_batch(W, rng, bs=32):
    idx = rng.choice(n_samples, bs, replace=False)
    Ai = A_opt[idx]
    Bi = B_opt[idx]
    return Ai.T @ (Ai @ W - Bi) / bs

def ema_smooth(data, alpha=0.95):
    """Exponential moving average smoothing"""
    result = np.zeros_like(data)
    result[0] = data[0]
    for i in range(1, len(data)):
        result[i] = alpha * result[i-1] + (1 - alpha) * data[i]
    return result

ns_iters_list = [1, 2, 5, 10]
colors_right = ['#e74c3c', '#f39c12', '#2ecc71', '#3498db']
n_steps = 3000
record_every = 3
mu_mom = 0.95
alpha_lr = 0.02

for ns_it, color in zip(ns_iters_list, colors_right):
    all_losses = []
    for seed in range(5):
        rng_opt = np.random.RandomState(seed)
        W = np.zeros((d_opt, p_opt))
        buf = np.zeros((d_opt, p_opt))
        losses = []
        for t in range(n_steps):
            if t % record_every == 0:
                losses.append(loss(W) - f_star)
            G = stoch_grad_batch(W, rng_opt)
            buf = mu_mom * buf + G
            G_tilde = mu_mom * buf + G
            O = newton_schulz(G_tilde, ns_it)
            W = W - alpha_lr * O
        all_losses.append(losses)
    median = np.median(all_losses, axis=0)
    smoothed = ema_smooth(np.maximum(median, 1e-16), alpha=0.97)
    steps = np.arange(0, n_steps, record_every)
    ax.semilogy(steps, smoothed, '-', color=color, linewidth=2.5,
                label=f'NS iters = {ns_it}', alpha=0.85)

# Exact SVD
all_losses_exact = []
for seed in range(5):
    rng_opt = np.random.RandomState(seed)
    W = np.zeros((d_opt, p_opt))
    buf = np.zeros((d_opt, p_opt))
    losses = []
    for t in range(n_steps):
        if t % record_every == 0:
            losses.append(loss(W) - f_star)
        G = stoch_grad_batch(W, rng_opt)
        buf = mu_mom * buf + G
        G_tilde = mu_mom * buf + G
        U_, _, Vt_ = np.linalg.svd(G_tilde, full_matrices=False)
        O = U_ @ Vt_
        W = W - alpha_lr * O
    all_losses_exact.append(losses)
median_exact = np.median(all_losses_exact, axis=0)
smoothed_exact = ema_smooth(np.maximum(median_exact, 1e-16), alpha=0.97)
ax.semilogy(steps, smoothed_exact, 'k--', linewidth=2,
            label='Exact SVD ($\\infty$ iters)', alpha=0.6)

ax.set_xlabel('Steps', fontsize=14)
ax.set_ylabel(r'$f(W^k) - f^*$', fontsize=13)
ax.set_title('Muon: effect of Newton-Schulz iterations', fontsize=14)
ax.legend(fontsize=11, loc='upper right')
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('/root/mipt25_work/files/exp_newton_schulz.pdf', bbox_inches='tight', dpi=150)
plt.savefig('/root/mipt25_work/files/exp_newton_schulz.png', bbox_inches='tight', dpi=150)
print("Done!")
