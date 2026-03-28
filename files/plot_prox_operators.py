import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Style
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['axes.linewidth'] = 1.2
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amssymb}'

y = np.linspace(-4, 4, 500)
lam = 1.0

# === Plot 1: Soft-thresholding (L1 prox) ===
def soft_threshold(y, lam):
    return np.sign(y) * np.maximum(np.abs(y) - lam, 0)

fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.0))

ax.plot(y, y, color='#BBBBBB', linewidth=1.2, linestyle='--', label=r'$\mathrm{Id}(y) = y$', zorder=1)
ax.plot(y, soft_threshold(y, lam), color='#E63946', linewidth=2.8, label=r'$\mathrm{prox}_{\lambda\|\cdot\|_1}(y),\ \lambda=1$', zorder=3)

ax.axhline(0, color='black', linewidth=0.5, zorder=0)
ax.axvline(0, color='black', linewidth=0.5, zorder=0)

# Mark lambda thresholds
ax.axvline(lam, color='#457B9D', linewidth=1.0, linestyle=':', alpha=0.7)
ax.axvline(-lam, color='#457B9D', linewidth=1.0, linestyle=':', alpha=0.7)
ax.annotate(r'$\lambda$', xy=(lam, 0), xytext=(lam + 0.15, -0.6), fontsize=13, color='#457B9D')
ax.annotate(r'$-\lambda$', xy=(-lam, 0), xytext=(-lam - 0.6, -0.6), fontsize=13, color='#457B9D')

# Dead zone shading
ax.axvspan(-lam, lam, alpha=0.08, color='#457B9D', zorder=0)

ax.set_xlabel(r'$y$', fontsize=15)
ax.set_ylabel(r'$\mathrm{prox}_r(y)$', fontsize=15)
ax.set_xlim(-4, 4)
ax.set_ylim(-3.5, 3.5)
ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
ax.set_aspect('equal')
ax.grid(True, alpha=0.2)

fig.tight_layout()
fig.savefig('/Users/bratishka/Yandex.Disk.localized/Base/Git/mipt25/files/prox_l1.pdf', bbox_inches='tight')
plt.close()

# === Plot 2: L2 prox (shrinkage) ===
fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.0))

lambdas = [0.5, 1.0, 2.0]
colors = ['#2A9D8F', '#E76F51', '#264653']

ax.plot(y, y, color='#BBBBBB', linewidth=1.2, linestyle='--', label=r'$\mathrm{Id}(y) = y$', zorder=1)

for lam_val, col in zip(lambdas, colors):
    prox_l2 = y / (1 + lam_val)
    ax.plot(y, prox_l2, color=col, linewidth=2.5,
            label=rf'$\lambda = {lam_val}$', zorder=3)

ax.axhline(0, color='black', linewidth=0.5, zorder=0)
ax.axvline(0, color='black', linewidth=0.5, zorder=0)

ax.set_xlabel(r'$y$', fontsize=15)
ax.set_ylabel(r'$\mathrm{prox}_r(y) = \dfrac{y}{1+\lambda}$', fontsize=15)
ax.set_xlim(-4, 4)
ax.set_ylim(-3.5, 3.5)
ax.legend(loc='upper left', fontsize=11, framealpha=0.9,
          title=r'$r(x)=\frac{\lambda}{2}\|x\|_2^2$')
ax.set_aspect('equal')
ax.grid(True, alpha=0.2)

fig.tight_layout()
fig.savefig('/Users/bratishka/Yandex.Disk.localized/Base/Git/mipt25/files/prox_l2.pdf', bbox_inches='tight')
plt.close()

# === Plot 3: Projection onto interval [-a, a] ===
a = 1.5

def proj_interval(y, a):
    return np.clip(y, -a, a)

fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.0))

ax.plot(y, y, color='#BBBBBB', linewidth=1.2, linestyle='--', label=r'$\mathrm{Id}(y) = y$', zorder=1)
ax.plot(y, proj_interval(y, a), color='#6A4C93', linewidth=2.8,
        label=rf'$\mathrm{{proj}}_{{[-a,a]}}(y),\ a={a}$', zorder=3)

ax.axhline(0, color='black', linewidth=0.5, zorder=0)
ax.axvline(0, color='black', linewidth=0.5, zorder=0)

# Mark interval boundaries
ax.axhline(a, color='#6A4C93', linewidth=1.0, linestyle=':', alpha=0.5)
ax.axhline(-a, color='#6A4C93', linewidth=1.0, linestyle=':', alpha=0.5)
ax.annotate(r'$a$', xy=(0, a), xytext=(-0.7, a + 0.15), fontsize=13, color='#6A4C93')
ax.annotate(r'$-a$', xy=(0, -a), xytext=(-0.9, -a - 0.4), fontsize=13, color='#6A4C93')

ax.set_xlabel(r'$y$', fontsize=15)
ax.set_ylabel(r'$\mathrm{prox}_{\mathbb{I}_S}(y) = \mathrm{proj}_S(y)$', fontsize=15)
ax.set_xlim(-4, 4)
ax.set_ylim(-3.5, 3.5)
ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
ax.set_aspect('equal')
ax.grid(True, alpha=0.2)

fig.tight_layout()
fig.savefig('/Users/bratishka/Yandex.Disk.localized/Base/Git/mipt25/files/prox_proj.pdf', bbox_inches='tight')
plt.close()

print("All plots saved successfully.")
