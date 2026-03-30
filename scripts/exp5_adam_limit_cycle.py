"""
Experiment 5: Adam does NOT converge on L(w) = w^2.

Reproduces Exercise 12, Ex. 2 (Adam non-convergence) from Orvieto's
"Nonconvex Optimization for Deep Learning" course.

Simplified Adam (beta1=0, eps=0) has a limit cycle at w^2 = eta^2/4
on the trivial quadratic L(w) = w^2. SGD and SGD+Momentum converge to 0.

Saves: exp_adam_limit_cycle.pdf, exp_adam_limit_cycle.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(42)

plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 16,
    'legend.fontsize': 11,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'mathtext.fontset': 'cm',
})

# ============================================================
# Parameters
# ============================================================
w0 = 1.0
eta = 0.1
n_iters = 2000

# L(w) = w^2, grad L(w) = 2w
def grad_L(w):
    return 2.0 * w

# ============================================================
# 1. SGD: w_{k+1} = w_k - eta * grad
# ============================================================
w_sgd = np.zeros(n_iters + 1)
w_sgd[0] = w0
for k in range(n_iters):
    w_sgd[k + 1] = w_sgd[k] - eta * grad_L(w_sgd[k])

# ============================================================
# 2. SGD + Momentum (heavy ball, mu=0.9)
# ============================================================
w_sgdm = np.zeros(n_iters + 1)
w_sgdm[0] = w0
vel = 0.0
for k in range(n_iters):
    vel = 0.9 * vel + grad_L(w_sgdm[k])
    w_sgdm[k + 1] = w_sgdm[k] - eta * vel

# ============================================================
# 3. Simplified Adam (beta1=0, eps=0, beta2=0.9)
#    Update (from Exercise 12):
#      v_k = beta2 * v_{k-1} + (1-beta2) * grad_k^2
#      w_{k+1} = w_k - eta * grad_k / sqrt(v_k)
#    Limit cycle: w^2 -> eta^2/4
# ============================================================
w_adam = np.zeros(n_iters + 1)
w_adam[0] = w0
v_adam = grad_L(w0)**2
for k in range(n_iters):
    g = grad_L(w_adam[k])
    v_adam = 0.9 * v_adam + 0.1 * g**2
    w_adam[k + 1] = w_adam[k] - eta * g / np.sqrt(v_adam)

# ============================================================
# Theoretical limit cycle: w^2 = eta^2 / 4
# ============================================================
limit_w_sq = eta**2 / 4   # 0.0025
limit_w = eta / 2         # 0.05

# Losses
L_sgd = w_sgd**2
L_sgdm = w_sgdm**2
L_adam = w_adam**2

iters = np.arange(n_iters + 1)

# Diagnostics
print(f"Adam final: w={w_adam[-1]:+.6f}, w^2={L_adam[-1]:.6f} "
      f"(theory eta^2/4={limit_w_sq:.6f})")
print(f"SGD final:  w={w_sgd[-1]:.2e}")
print(f"SGD+M final: w={w_sgdm[-1]:.2e}")

# ============================================================
# Plot
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle(r"Adam не сходится на $L(w) = w^2$",
             fontsize=18, fontweight='bold')

# ====== Left panel: w_k trajectory ======
ax1.plot(iters, w_sgd, color='black', linewidth=2.5,
         label='SGD', alpha=0.85)
ax1.plot(iters, w_sgdm, color='#e67e22', linewidth=2.5,
         label=r'SGD + Momentum ($\mu\!=\!0.9$)', alpha=0.85)
ax1.plot(iters, w_adam, color='#e74c3c', linewidth=0.8,
         label=r'Adam ($\beta_1\!=\!0,\;\beta_2\!=\!0.9$)', alpha=0.85)

ax1.axhline(y=limit_w, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)
ax1.axhline(y=-limit_w, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)
ax1.text(n_iters * 0.78, limit_w + 0.02, r'$\eta/2$',
         color='gray', fontsize=14, ha='center')
ax1.text(n_iters * 0.78, -limit_w - 0.04, r'$-\eta/2$',
         color='gray', fontsize=14, ha='center')

ax1.set_xlabel('Итерация')
ax1.set_ylabel(r'Значение параметра $w_k$')
ax1.legend(loc='upper right', framealpha=0.95, edgecolor='lightgray')
ax1.grid(True, alpha=0.25)
ax1.set_xlim(0, n_iters)

# Inset: zoom into a narrow window to see individual oscillations
axins = ax1.inset_axes([0.35, 0.05, 0.62, 0.42])
t_start, t_end = 1000, 1040
axins.plot(iters[t_start:t_end], w_adam[t_start:t_end], color='#e74c3c',
           linewidth=1.2, marker='o', markersize=3, alpha=0.85)
axins.axhline(y=limit_w, color='gray', linestyle='--', linewidth=1.2, alpha=0.6)
axins.axhline(y=-limit_w, color='gray', linestyle='--', linewidth=1.2, alpha=0.6)
axins.axhline(y=0, color='lightgray', linewidth=0.5)
axins.set_xlim(t_start, t_end)
axins.set_ylim(-limit_w * 1.8, limit_w * 1.8)
axins.set_ylabel(r'$w_k$', fontsize=11)
axins.set_xticks([1000, 1010, 1020, 1030, 1040])
axins.tick_params(labelsize=9)
axins.grid(True, alpha=0.2)
axins.set_title(r'Увеличено: $w_k$ осциллирует около $\pm\eta/2$', fontsize=10,
                pad=3)

# ====== Right panel: L(w_k) = w_k^2 ======
floor = 1e-320
ax2.semilogy(iters, L_sgd + floor, color='black', linewidth=2.5,
             label='SGD', alpha=0.85)
ax2.semilogy(iters, L_sgdm + floor, color='#e67e22', linewidth=2.5,
             label=r'SGD + Momentum', alpha=0.85)
ax2.semilogy(iters, L_adam + floor, color='#e74c3c', linewidth=1.5,
             label='Adam', alpha=0.85)

ax2.axhline(y=limit_w_sq, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.text(n_iters * 0.55, limit_w_sq * 5,
         r'$\eta^2/4 = %.4f$' % limit_w_sq,
         color='gray', fontsize=13, ha='center')

ax2.set_xlabel('Итерация')
ax2.set_ylabel(r'Значение функции $L(w_k) = w_k^2$')
ax2.legend(loc='center right', framealpha=0.95, edgecolor='lightgray')
ax2.grid(True, alpha=0.25)
ax2.set_xlim(0, n_iters)

fig.tight_layout()

# ============================================================
# Save
# ============================================================
out = '/root/mipt25_work/files'
fig.savefig(f'{out}/exp_adam_limit_cycle.pdf', bbox_inches='tight', dpi=150)
fig.savefig(f'{out}/exp_adam_limit_cycle.png', bbox_inches='tight', dpi=150)
plt.close()
print(f"\nSaved to {out}/exp_adam_limit_cycle.{{pdf,png}}")
