"""
Experiment 4: Muon geometry — what orthogonalization does to gradient directions.
Three panels:
1. Singular values of G vs UV^T (all ones for Muon)
2. Update energy per direction for SGD/Adam/Muon
3. Condition number of gradient vs Muon direction over training
All labels in Russian. Publication quality.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

np.random.seed(42)

plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 16,
    'legend.fontsize': 11,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

fig, axes = plt.subplots(1, 3, figsize=(17, 5.5),
                         gridspec_kw={'width_ratios': [1, 1, 1.2], 'wspace': 0.36})
ax1, ax2, ax3 = axes

# ============================================================
# Gradient matrix with very different singular values
# ============================================================
d1, d2 = 8, 6
rng = np.random.RandomState(42)
U = np.linalg.qr(rng.randn(d1, d1))[0][:, :d2]
V = np.linalg.qr(rng.randn(d2, d2))[0]
sigma = np.array([50, 20, 5, 2, 0.5, 0.1])

# ============================================================
# Panel 1: Singular values of G vs UV^T
# ============================================================
x = np.arange(len(sigma))
width = 0.35

ax1.bar(x - width/2, sigma, width, color='#e74c3c', alpha=0.8,
        label='$\\sigma_i(G)$', edgecolor='white')
ax1.bar(x + width/2, np.ones(len(sigma)), width, color='#2980b9', alpha=0.8,
        label='$\\sigma_i(UV^\\top) = 1$', edgecolor='white')

ax1.set_xlabel('Индекс $i$', fontsize=15)
ax1.set_ylabel('Сингулярное число', fontsize=15)
ax1.set_title('Сингулярные числа\nградиента и Muon', fontsize=16, pad=12)
ax1.legend(fontsize=11, loc='center right')
ax1.set_yscale('log')
ax1.set_xticks(x)
ax1.set_xticklabels([f'{i+1}' for i in x])
ax1.grid(True, alpha=0.3, axis='y')
ax1.tick_params(labelsize=12)
ax1.set_ylim(0.05, 200)

# Annotation: point from sigma_1 bar, text placed clearly to the right
ax1.annotate('разница в 500 раз', xy=(0.15, 50), xytext=(3.2, 120),
             fontsize=11, color='#c0392b', fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5),
             ha='center', va='center')

# ============================================================
# Panel 2: Update energy per direction for SGD/Adam/Muon
# ============================================================
sgd_update_norms = sigma
adam_update_norms = np.sqrt(sigma)
muon_update_norms = np.ones(len(sigma))

x = np.arange(len(sigma))
width = 0.25

ax2.bar(x - width, sgd_update_norms / sgd_update_norms.max(), width,
        color='#e67e22', alpha=0.8, label='SGD', edgecolor='white')
ax2.bar(x, adam_update_norms / adam_update_norms.max(), width,
        color='#27ae60', alpha=0.8, label='Adam (прибл.)', edgecolor='white')
ax2.bar(x + width, muon_update_norms, width,
        color='#2980b9', alpha=0.8, label='Muon', edgecolor='white')

ax2.set_xlabel('Индекс направления $i$', fontsize=15)
ax2.set_ylabel('Отн. величина обновления', fontsize=15)
ax2.set_title('Энергия обновления\nпо направлениям', fontsize=16, pad=12)
ax2.legend(fontsize=11, loc='upper right')
ax2.set_xticks(x)
ax2.set_xticklabels([f'{i+1}' for i in x])
ax2.grid(True, alpha=0.3, axis='y')
ax2.tick_params(labelsize=12)

# ============================================================
# Panel 3: Condition number over training
# ============================================================
m_opt, d_opt, p_opt = 100, 15, 10
rng = np.random.RandomState(42)
U_A = np.linalg.qr(rng.randn(m_opt, m_opt))[0][:, :d_opt]
V_A = np.linalg.qr(rng.randn(d_opt, d_opt))[0]
s_A = np.logspace(0, 3, d_opt)
A = U_A @ np.diag(s_A) @ V_A.T
B = rng.randn(m_opt, p_opt) * 0.1

def stoch_grad(W, rng_loc, bs=16):
    idx = rng_loc.choice(m_opt, bs, replace=False)
    return A[idx].T @ (A[idx] @ W - B[idx]) / bs

W = np.zeros((d_opt, p_opt))
buf = np.zeros((d_opt, p_opt))
cond_grad = []

n_track = 500
for t in range(n_track):
    G_t = stoch_grad(W, rng, bs=16)
    buf = 0.95 * buf + G_t
    G_tilde = 0.95 * buf + G_t

    s_g = np.linalg.svd(G_tilde, compute_uv=False)
    cond_g = s_g[0] / (s_g[-1] + 1e-12)
    cond_grad.append(cond_g)

    # Newton-Schulz polar decomposition
    a, b, c = (3.4445, -4.7750, 2.0315)
    X_ns = G_tilde / (np.linalg.norm(G_tilde, 'fro') + 1e-12)
    X_ns = X_ns / (max(G_tilde.shape)**0.5)
    for _ in range(5):
        A_ = X_ns @ X_ns.T
        B_ = b * A_ + c * A_ @ A_
        X_ns = a * X_ns + B_ @ X_ns
    W = W - 0.02 * X_ns

cond_grad_smooth = uniform_filter1d(np.log10(cond_grad), size=20)

ax3.plot(range(n_track), 10**cond_grad_smooth, '-', color='#e74c3c', linewidth=2.5,
         label=r'$\mathrm{cond}(G)$ — градиент', alpha=0.85)
ax3.axhline(y=1, color='#2980b9', linewidth=2.5, linestyle='-',
            label=r'$\mathrm{cond}(UV^\top) = 1$ — Muon', alpha=0.85)

ax3.set_xlabel('Шаг обучения', fontsize=15)
ax3.set_ylabel('Число обусловленности', fontsize=15)
ax3.set_title('Число обусловленности', fontsize=16, pad=12)
ax3.legend(fontsize=11, loc='upper right')
ax3.set_yscale('log')
ax3.grid(True, alpha=0.3)
ax3.tick_params(labelsize=12)

# Shaded region
ax3.fill_between(range(n_track), 1, 10**cond_grad_smooth, alpha=0.08, color='#e74c3c')

# Place annotation text in the middle of the shaded region
mid_log = np.mean(cond_grad_smooth) * 0.5
ax3.text(250, 10**mid_log, 'Muon устраняет\nэту анизотропию',
         fontsize=12, ha='center', color='#8e44ad', alpha=0.85, fontweight='bold')

ymax_log = np.max(cond_grad_smooth) + 0.5
ax3.set_ylim(0.5, 10**ymax_log)

fig.subplots_adjust(left=0.04, right=0.98, bottom=0.13, top=0.88, wspace=0.36)
plt.savefig('/root/mipt25_work/files/exp_muon_geometry.pdf', bbox_inches='tight', dpi=150)
plt.savefig('/root/mipt25_work/files/exp_muon_geometry.png', bbox_inches='tight', dpi=150)
print("Saved exp_muon_geometry.pdf and .png")
