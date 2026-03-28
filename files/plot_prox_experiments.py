"""
Numerical experiments for Lecture 18: Proximal Gradient Methods.
Generates PDF plots for beamer slides.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import cvxpy as cp
import warnings
warnings.filterwarnings('ignore')

# ── Style ────────────────────────────────────────────────────────────────
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 13
rcParams['axes.linewidth'] = 1.2
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = (
    r'\usepackage{amsmath}\usepackage{amssymb}'
    r'\usepackage[T2A]{fontenc}\usepackage[utf8]{inputenc}\usepackage[russian]{babel}'
)

OUT = '/Users/bratishka/Yandex.Disk.localized/Base/Git/mipt25/files'

# Colors
C_SUBGRAD = '#E63946'
C_PROX    = '#2A9D8F'
C_FISTA   = '#264653'
C_THEORY  = '#999999'
C_STAR    = '#457B9D'

np.random.seed(42)


# ── Utilities ────────────────────────────────────────────────────────────
def soft_threshold(x, kappa):
    return np.sign(x) * np.maximum(np.abs(x) - kappa, 0.0)


def svd_soft_threshold(X, kappa):
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    s_new = np.maximum(s - kappa, 0.0)
    return (U * s_new[None, :]) @ Vt


def generate_A(m, n, mu_H, L_H):
    """Generate A such that (1/m)*A^T A has eigenvalues in [mu_H, L_H]."""
    eigs = np.linspace(mu_H, L_H, n)
    Q, _ = np.linalg.qr(np.random.randn(n, n))
    # A = random_m_x_n @ Q @ diag(sqrt(eigs * m)) @ Q^T → (1/m)A^T A = Q diag(eigs) Q^T
    # Simpler: create A with prescribed singular values
    U, _ = np.linalg.qr(np.random.randn(m, n), mode='reduced')  # m x n orthonormal
    svals = np.sqrt(eigs * m)
    A = U @ np.diag(svals) @ Q.T
    return A


# ── Algorithms ───────────────────────────────────────────────────────────
def subgradient_method(x0, subgrad_fn, phi_fn, step_fn, n_iter):
    """Subgradient with averaged iterate."""
    x = x0.copy()
    x_sum = x0.copy()
    phi_avg = [phi_fn(x0)]
    x_norms = [0.0]
    G_max = 0.0
    for k in range(1, n_iter + 1):
        g = subgrad_fn(x)
        G_max = max(G_max, np.linalg.norm(g))
        alpha_k = step_fn(k, G_max)
        x = x - alpha_k * g
        x_sum = x_sum + x
        x_avg = x_sum / (k + 1)
        phi_avg.append(phi_fn(x_avg))
    return np.array(phi_avg), G_max


def subgradient_method_xgap(x0, subgrad_fn, x_star, step_fn, n_iter):
    """Subgradient tracking ||x_avg - x*||^2."""
    x = x0.copy()
    x_sum = x0.copy()
    gaps = [np.linalg.norm(x0 - x_star)**2]
    G_max = 0.0
    for k in range(1, n_iter + 1):
        g = subgrad_fn(x)
        G_max = max(G_max, np.linalg.norm(g))
        alpha_k = step_fn(k, G_max)
        x = x - alpha_k * g
        x_sum = x_sum + x
        x_avg = x_sum / (k + 1)
        gaps.append(np.linalg.norm(x_avg - x_star)**2)
    return np.array(gaps), G_max


def proximal_gradient(x0, grad_f, prox_r, alpha, phi_fn, n_iter):
    """ISTA."""
    x = x0.copy()
    vals = [phi_fn(x0)]
    for _ in range(n_iter):
        x = prox_r(x - alpha * grad_f(x), alpha)
        vals.append(phi_fn(x))
    return np.array(vals), x


def proximal_gradient_xgap(x0, grad_f, prox_r, alpha, x_star, n_iter):
    """ISTA tracking ||x_k - x*||^2."""
    x = x0.copy()
    gaps = [np.linalg.norm(x0 - x_star)**2]
    for _ in range(n_iter):
        x = prox_r(x - alpha * grad_f(x), alpha)
        gaps.append(np.linalg.norm(x - x_star)**2)
    return np.array(gaps), x


def fista(x0, grad_f, prox_r, alpha, phi_fn, n_iter):
    """FISTA (accelerated proximal gradient)."""
    x = x0.copy()
    y = x0.copy()
    t = 1.0
    vals = [phi_fn(x0)]
    for _ in range(n_iter):
        x_new = prox_r(y - alpha * grad_f(y), alpha)
        t_new = (1 + np.sqrt(1 + 4 * t**2)) / 2
        y = x_new + ((t - 1) / t_new) * (x_new - x)
        x, t = x_new, t_new
        vals.append(phi_fn(x))
    return np.array(vals), x


# ══════════════════════════════════════════════════════════════════════════
# Plot 1: Convex LASSO — theoretical rates
# ══════════════════════════════════════════════════════════════════════════
def plot_convex_theory():
    print("Plot 1: Convex LASSO...")
    m, n, lam, L_H = 200, 100, 0.5, 10.0
    mu_H = 0.01  # nearly zero for convex

    A = generate_A(m, n, mu_H, L_H)
    w_true = np.random.randn(n)
    w_true[np.random.choice(n, n // 2, replace=False)] = 0
    b = A @ w_true + 0.1 * np.random.randn(m)

    L = L_H  # smoothness of f

    # Optimal via CVXPY
    w_var = cp.Variable(n)
    prob = cp.Problem(cp.Minimize(0.5 / m * cp.sum_squares(A @ w_var - b) + lam * cp.norm1(w_var)))
    prob.solve(solver=cp.SCS, verbose=False)
    x_star, f_star = w_var.value, prob.value

    # Functions
    grad_f = lambda x: (1.0 / m) * A.T @ (A @ x - b)
    phi = lambda x: 0.5 / m * np.linalg.norm(A @ x - b)**2 + lam * np.linalg.norm(x, 1)
    subgrad = lambda x: grad_f(x) + lam * np.sign(x)
    prox_l1 = lambda x, alpha: soft_threshold(x, lam * alpha)

    x0 = 0.5 * np.random.randn(n)
    R = np.linalg.norm(x0 - x_star)
    N = 5000

    # Subgradient: step = R / (G * sqrt(k)), G estimated online
    step_sub = lambda k, G: R / (max(G, 1e-8) * np.sqrt(k))
    phi_sub, G = subgradient_method(x0, subgrad, phi, step_sub, N)

    # Proximal
    phi_prox, _ = proximal_gradient(x0, grad_f, prox_l1, 1.0 / L, phi, N)

    # FISTA
    phi_fista, _ = fista(x0, grad_f, prox_l1, 1.0 / L, phi, N)

    # Gaps
    gap_sub = np.maximum(phi_sub - f_star, 1e-16)
    gap_prox = np.maximum(phi_prox - f_star, 1e-16)
    gap_fista = np.maximum(phi_fista - f_star, 1e-16)

    k = np.arange(1, N + 2)

    # Theoretical bounds
    bound_sub = G * R / np.sqrt(k)
    bound_prox = L * R**2 / (2 * k)
    bound_fista = 2 * L * R**2 / (k)**2

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.loglog(k, gap_sub, color=C_SUBGRAD, lw=2.2, label=r'Субградиентный ($\alpha_k = \frac{R}{G\sqrt{k}}$)')
    ax.loglog(k, gap_prox, color=C_PROX, lw=2.2, label=r'Проксимальный ($\alpha = \frac{1}{L}$)')
    ax.loglog(k, gap_fista, color=C_FISTA, lw=2.2, label=r'FISTA ($\alpha = \frac{1}{L}$)')

    ax.loglog(k, bound_sub, '--', color=C_SUBGRAD, alpha=0.5, lw=1.5, label=r'$\frac{GR}{\sqrt{k}}$')
    ax.loglog(k, bound_prox, '--', color=C_PROX, alpha=0.5, lw=1.5, label=r'$\frac{LR^2}{2k}$')
    ax.loglog(k, bound_fista, '--', color=C_FISTA, alpha=0.5, lw=1.5, label=r'$\frac{2LR^2}{k^2}$')

    # Starting point
    # No starting point annotation

    ax.set_xlabel(r'Итерации, $k$', fontsize=14)
    ax.set_ylabel(r'$\varphi(x_k) - \varphi^*$', fontsize=14)
    ax.set_title(r'LASSO: выпуклый случай ($\mu \approx 0$). $m=%d$, $n=%d$, $\lambda=%.1f$, $L=%d$' % (m, n, lam, int(L)), fontsize=13)
    ax.legend(fontsize=9, ncol=2, loc='lower right', framealpha=0.9)
    ax.grid(True, alpha=0.2, which='both')
    ax.set_xlim(1, N)

    fig.tight_layout()
    fig.savefig(f'{OUT}/prox_convex_theory.pdf', bbox_inches='tight')
    plt.close()
    print(f"  R={R:.2f}, G={G:.2f}, L={L:.1f}")


# ══════════════════════════════════════════════════════════════════════════
# Plot 2: Strongly convex LASSO — theoretical rates
# ══════════════════════════════════════════════════════════════════════════
def plot_sc_theory():
    print("Plot 2: Strongly convex LASSO...")
    m, n, lam, mu, L = 200, 100, 0.5, 1.0, 10.0

    # Eigenvalues of (1/m)*A^T A in [0, L-mu], so f(x) = LS + (mu/2)||x||^2 has spectrum [mu, L]
    A = generate_A(m, n, 0.01, L - mu)
    w_true = np.random.randn(n)
    w_true[np.random.choice(n, n // 2, replace=False)] = 0
    b = A @ w_true + 0.1 * np.random.randn(m)

    # f(x) = (1/2m)||Ax-b||^2 + (mu/2)||x||^2   (smooth, mu-strongly convex, L-smooth)
    grad_f = lambda x: (1.0 / m) * A.T @ (A @ x - b) + mu * x
    f = lambda x: 0.5 / m * np.linalg.norm(A @ x - b)**2 + 0.5 * mu * np.linalg.norm(x)**2
    phi = lambda x: f(x) + lam * np.linalg.norm(x, 1)
    subgrad = lambda x: grad_f(x) + lam * np.sign(x)
    prox_l1 = lambda x, alpha: soft_threshold(x, lam * alpha)

    # Optimal
    w_var = cp.Variable(n)
    prob = cp.Problem(cp.Minimize(
        0.5 / m * cp.sum_squares(A @ w_var - b) + 0.5 * mu * cp.sum_squares(w_var) + lam * cp.norm1(w_var)
    ))
    prob.solve(solver=cp.SCS, verbose=False, max_iters=100000, eps=1e-12)
    # Polish x* by running many proximal iterations from CVXPY solution
    x_star_polished = w_var.value.copy()
    for _ in range(50000):
        x_star_polished = prox_l1(x_star_polished - (1.0/L) * grad_f(x_star_polished), 1.0/L)
    x_star = x_star_polished

    x0 = 0.5 * np.random.randn(n)
    R2 = np.linalg.norm(x0 - x_star)**2
    N = 500

    # Subgradient: step = 2/(mu*(k+1))
    step_sub = lambda k, G: 2.0 / (mu * (k + 1))
    xgap_sub, G = subgradient_method_xgap(x0, subgrad, x_star, step_sub, N)

    # Proximal
    alpha = 1.0 / L
    xgap_prox, _ = proximal_gradient_xgap(x0, grad_f, prox_l1, alpha, x_star, N)

    k = np.arange(0, N + 1)

    # Theoretical bounds
    bound_sub = 2 * G**2 / (mu**2 * np.maximum(k, 1))
    bound_prox = (1 - mu / L)**k * R2

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.semilogy(k, xgap_sub, color=C_SUBGRAD, lw=2.2,
                label=r'Субградиентный ($\alpha_k = \frac{2}{\mu(k+1)}$)')
    ax.semilogy(k, xgap_prox, color=C_PROX, lw=2.2,
                label=r'Проксимальный ($\alpha = \frac{1}{L}$)')

    ax.semilogy(k[1:], bound_sub[1:], '--', color=C_SUBGRAD, alpha=0.5, lw=1.5,
                label=r'$\frac{2G^2}{\mu^2 k}$')
    ax.semilogy(k, bound_prox, '--', color=C_PROX, alpha=0.5, lw=1.5,
                label=r'$(1-\frac{\mu}{L})^k \|x_0-x^*\|^2$')

    # Starting point
    # No starting point annotation

    ax.set_xlabel(r'Итерации, $k$', fontsize=14)
    ax.set_ylabel(r'$\|x_k - x^*\|^2$', fontsize=14)
    ax.set_title(r'LASSO: сильно выпуклый случай. $\mu=%g$, $L=%d$, $\kappa = %d$' % (mu, int(L), int(L/mu)), fontsize=13)
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
    ax.grid(True, alpha=0.2, which='both')
    ax.set_xlim(0, N)

    fig.tight_layout()
    fig.savefig(f'{OUT}/prox_sc_theory.pdf', bbox_inches='tight')
    plt.close()
    print(f"  R^2={R2:.2f}, G={G:.2f}, mu/L={mu/L:.2f}")


# ══════════════════════════════════════════════════════════════════════════
# Plot 3: Sparsity emergence
# ══════════════════════════════════════════════════════════════════════════
def plot_sparsity():
    print("Plot 3: Sparsity...")
    m, n, lam, L_H = 200, 30, 1.5, 10.0

    A = generate_A(m, n, 0.01, L_H)
    w_true = np.zeros(n)
    w_true[:10] = np.random.randn(10) * 2
    b = A @ w_true + 0.1 * np.random.randn(m)

    L = L_H
    grad_f = lambda x: (1.0 / m) * A.T @ (A @ x - b)
    phi = lambda x: 0.5 / m * np.linalg.norm(A @ x - b)**2 + lam * np.linalg.norm(x, 1)
    subgrad = lambda x: grad_f(x) + lam * np.sign(x)
    prox_l1 = lambda x, alpha: soft_threshold(x, lam * alpha)

    # Optimal — polish with proximal iterations for exact zeros
    w_var = cp.Variable(n)
    prob = cp.Problem(cp.Minimize(0.5 / m * cp.sum_squares(A @ w_var - b) + lam * cp.norm1(w_var)))
    prob.solve(solver=cp.SCS, verbose=False)
    x_star = w_var.value.copy()
    for _ in range(20000):
        x_star = prox_l1(x_star - (1.0 / L) * grad_f(x_star), 1.0 / L)

    x0 = 0.3 * np.random.randn(n)
    R = np.linalg.norm(x0 - x_star)
    N = 3000

    # Subgradient
    x_sub = x0.copy()
    x_sub_sum = x0.copy()
    G_max = 0.0
    for k in range(1, N + 1):
        g = subgrad(x_sub)
        G_max = max(G_max, np.linalg.norm(g))
        alpha_k = R / (max(G_max, 1e-8) * np.sqrt(k))
        x_sub = x_sub - alpha_k * g
        x_sub_sum += x_sub
    x_sub_final = x_sub_sum / (N + 1)

    # Proximal
    x_prox = x0.copy()
    for _ in range(N):
        x_prox = prox_l1(x_prox - (1.0 / L) * grad_f(x_prox), 1.0 / L)

    thr = 1e-10
    n_zero_star = np.sum(np.abs(x_star) < thr)
    n_zero_sub = np.sum(np.abs(x_sub_final) < thr)
    n_zero_prox = np.sum(np.abs(x_prox) < thr)

    # Identify "zero" components of x*
    zero_mask = np.abs(x_star) < thr
    zero_idx = np.where(zero_mask)[0]
    nonzero_idx = np.where(~zero_mask)[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2),
                                    gridspec_kw={'width_ratios': [1, 1.2]})

    # Left: full solution on linear scale
    idx = np.arange(n)
    w = 0.3
    ax1.bar(idx - w/2, x_sub_final, w, color=C_SUBGRAD, alpha=0.8, label=r'Субградиентный')
    ax1.bar(idx + w/2, x_prox, w, color=C_PROX, alpha=0.8, label=r'Проксимальный')
    ax1.scatter(idx, x_star, s=18, color=C_STAR, zorder=5, marker='o', label=r'$x^*$', edgecolors='white', linewidths=0.3)
    ax1.set_xlabel(r'Компонента $i$', fontsize=11)
    ax1.set_ylabel(r'$x_i$', fontsize=12)
    ax1.set_title(r'Полное решение', fontsize=12)
    ax1.legend(fontsize=8, loc='upper right', framealpha=0.9)
    ax1.grid(True, alpha=0.15)

    # Right: zoom on "zero" components
    zidx = np.arange(len(zero_idx))
    ax2.bar(zidx - w/2, x_sub_final[zero_idx], w, color=C_SUBGRAD, alpha=0.8,
            label=r'Субград. ($|x_i| \approx 10^{-2}$)')
    ax2.bar(zidx + w/2, x_prox[zero_idx], w, color=C_PROX, alpha=0.8,
            label=r'Прокс. ($x_i = 0$ точно)')
    ax2.axhline(0, color='black', lw=0.5)
    ax2.set_xlabel(r'Компоненты, где $x^*_i = 0$', fontsize=11)
    ax2.set_ylabel(r'$x_i$', fontsize=12)
    ax2.set_title(r'Увеличение: «нулевые» компоненты', fontsize=12)
    ax2.legend(fontsize=8.5, loc='upper right', framealpha=0.9)
    ax2.grid(True, alpha=0.15)
    ax2.set_xticks(zidx[::3])
    ax2.set_xticklabels(zero_idx[::3])

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle(r'LASSO: разреженность. $n=%d$, $\lambda=%.1f$. Нулей: субград. %d/%d, прокс. %d/%d'
                 % (n, lam, n_zero_sub, n, n_zero_prox, n), fontsize=11)
    fig.savefig(f'{OUT}/prox_sparsity.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Zeros: x*={n_zero_star}, subgrad={n_zero_sub}, prox={n_zero_prox}")


# ══════════════════════════════════════════════════════════════════════════
# Plot 4: Binary logistic regression with L1
# ══════════════════════════════════════════════════════════════════════════
def plot_logistic():
    print("Plot 4: Logistic regression...")
    from sklearn.datasets import make_classification
    from sklearn.preprocessing import StandardScaler

    m, n_feat, lam = 300, 50, 0.1
    X, y01 = make_classification(n_samples=m, n_features=n_feat, n_informative=40,
                                  n_redundant=5, random_state=42)
    X = StandardScaler().fit_transform(X)
    y = 2 * y01 - 1  # {-1, +1}

    # L for logistic: L = (1/(4m)) * ||X^T X||_op
    L = np.linalg.norm(X.T @ X, ord=2) / (4 * m)

    def f(w):
        z = y * (X @ w)
        return np.mean(np.logaddexp(0, -z))

    def grad_f(w):
        z = y * (X @ w)
        s = -y / (1 + np.exp(z))  # -y * sigmoid(-z)
        return (1.0 / m) * X.T @ s

    phi = lambda w: f(w) + lam * np.linalg.norm(w, 1)
    subgrad = lambda w: grad_f(w) + lam * np.sign(w)
    prox_l1 = lambda w, alpha: soft_threshold(w, lam * alpha)

    # Optimal
    w_var = cp.Variable(n_feat)
    z_cp = cp.multiply(y, X @ w_var)
    prob = cp.Problem(cp.Minimize(cp.sum(cp.logistic(-z_cp)) / m + lam * cp.norm1(w_var)))
    prob.solve(solver=cp.SCS, verbose=False)
    x_star, f_star = w_var.value, prob.value

    x0 = np.zeros(n_feat)
    R = np.linalg.norm(x0 - x_star)
    N = 300

    step_sub = lambda k, G: R / (max(G, 1e-8) * np.sqrt(k))
    phi_sub, G = subgradient_method(x0, subgrad, phi, step_sub, N)
    phi_prox, _ = proximal_gradient(x0, grad_f, prox_l1, 1.0 / L, phi, N)
    phi_fista, _ = fista(x0, grad_f, prox_l1, 1.0 / L, phi, N)

    gap_sub = np.maximum(phi_sub - f_star, 1e-16)
    gap_prox = np.maximum(phi_prox - f_star, 1e-16)
    gap_fista = np.maximum(phi_fista - f_star, 1e-16)

    k = np.arange(0, N + 1)

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.semilogy(k, gap_sub, color=C_SUBGRAD, lw=2.2, label=r'Субградиентный')
    ax.semilogy(k, gap_prox, color=C_PROX, lw=2.2, label=r'Проксимальный (ISTA)')
    ax.semilogy(k, gap_fista, color=C_FISTA, lw=2.2, label=r'FISTA')

    # No starting point annotation

    ax.set_xlabel(r'Итерации, $k$', fontsize=14)
    ax.set_ylabel(r'$\varphi(x_k) - \varphi^*$', fontsize=14)
    ax.set_title(r'Логистическая регрессия с $\ell_1$. $m=%d$, $n=%d$, $\lambda=%.1f$' % (m, n_feat, lam), fontsize=13)
    ax.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax.grid(True, alpha=0.2, which='both')
    ax.set_xlim(0, N)

    fig.tight_layout()
    fig.savefig(f'{OUT}/prox_logistic.pdf', bbox_inches='tight')
    plt.close()
    print(f"  L={L:.4f}, R={R:.2f}, G={G:.2f}")


# ══════════════════════════════════════════════════════════════════════════
# Plot 5: Matrix completion with nuclear norm
# ══════════════════════════════════════════════════════════════════════════
def plot_matrix_completion():
    print("Plot 5: Matrix completion...")
    n1, n2, r_true = 50, 50, 5
    lam = 0.5
    p_observe = 0.5

    # Ground truth low-rank matrix
    U_true = np.random.randn(n1, r_true)
    V_true = np.random.randn(n2, r_true)
    M_true = U_true @ V_true.T

    # Observation mask
    Omega = np.random.rand(n1, n2) < p_observe
    M_obs = M_true * Omega

    # f(X) = 0.5 * ||P_Omega(X) - P_Omega(M)||_F^2
    # grad f(X) = P_Omega(X - M) = (X - M) * Omega
    # L = 1
    L_mc = 1.0

    def grad_f(X):
        return (X - M_true) * Omega

    def phi(X):
        residual = (X - M_true) * Omega
        return 0.5 * np.linalg.norm(residual, 'fro')**2 + lam * np.linalg.norm(X, 'nuc')

    alpha = 1.0 / L_mc

    # Optimal via CVXPY
    X_var = cp.Variable((n1, n2))
    residual_cp = cp.multiply(Omega.astype(float), X_var - M_true)
    prob = cp.Problem(cp.Minimize(0.5 * cp.sum_squares(residual_cp) + lam * cp.normNuc(X_var)))
    prob.solve(solver=cp.SCS, verbose=False)
    f_star = prob.value
    X_star = X_var.value

    # Proximal gradient for nuclear norm
    X = np.zeros((n1, n2))
    N = 500
    phi_vals = [phi(X)]
    ranks = [np.linalg.matrix_rank(X, tol=1e-4)]
    rel_errors = [np.linalg.norm(X - M_true, 'fro') / np.linalg.norm(M_true, 'fro')]

    for _ in range(N):
        X = svd_soft_threshold(X - alpha * grad_f(X), lam * alpha)
        phi_vals.append(phi(X))
        ranks.append(np.linalg.matrix_rank(X, tol=1e-4))
        rel_errors.append(np.linalg.norm(X - M_true, 'fro') / np.linalg.norm(M_true, 'fro'))

    phi_vals = np.array(phi_vals)
    gap = np.maximum(phi_vals - f_star, 1e-16)

    R2 = np.linalg.norm(np.zeros((n1, n2)) - X_star, 'fro')**2
    k = np.arange(0, N + 1)
    bound = L_mc * R2 / (2 * np.maximum(k, 1))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2))

    # Left: convergence
    ax1.semilogy(k, gap, color=C_PROX, lw=2.2, label=r'Проксимальный метод')
    ax1.semilogy(k[1:], bound[1:], '--', color=C_PROX, alpha=0.5, lw=1.5, label=r'$\frac{LR^2}{2k}$')
    # No starting point annotation
    ax1.set_xlabel(r'Итерации, $k$', fontsize=12)
    ax1.set_ylabel(r'$\varphi(X_k) - \varphi^*$', fontsize=12)
    ax1.set_title(r'Сходимость', fontsize=12)
    ax1.legend(fontsize=9, framealpha=0.9)
    ax1.grid(True, alpha=0.2, which='both')

    # Right: relative error + rank
    color_err = C_FISTA
    ax2.plot(k, rel_errors, color=color_err, lw=2.2, label=r'$\|X_k - M\|_F / \|M\|_F$')
    ax2.set_xlabel(r'Итерации, $k$', fontsize=12)
    ax2.set_ylabel(r'Отн. ошибка восстановления', fontsize=11, color=color_err)
    ax2.tick_params(axis='y', labelcolor=color_err)
    ax2.set_title(r'Восстановление матрицы', fontsize=12)
    ax2.grid(True, alpha=0.2)

    ax2r = ax2.twinx()
    ax2r.plot(k[1:], ranks[1:], color=C_SUBGRAD, lw=1.5, ls='--', alpha=0.7, label=r'$\mathrm{rank}(X_k)$')
    ax2r.axhline(r_true, color='gray', ls=':', lw=1, alpha=0.5)
    ax2r.set_ylabel(r'Ранг $X_k$', fontsize=11, color=C_SUBGRAD)
    ax2r.tick_params(axis='y', labelcolor=C_SUBGRAD)
    ax2r.annotate(r'Истинный ранг $r=%d$' % r_true, xy=(N * 0.6, r_true),
                  xytext=(N * 0.55, r_true + 8), fontsize=9, color='gray',
                  arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle(r'Восстановление матриц: ядерная норма. $%d \times %d$, $r=%d$, $\lambda=%.1f$, наблюдаем $%.0f\%%$'
                 % (n1, n2, r_true, lam, p_observe * 100), fontsize=11)
    fig.savefig(f'{OUT}/prox_matrix_completion.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Final rank: {ranks[-1]}, rel error: {rel_errors[-1]:.4f}")


# ══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    plot_convex_theory()
    plot_sc_theory()
    plot_sparsity()
    plot_logistic()
    plot_matrix_completion()
    print("\nAll plots saved to", OUT)
