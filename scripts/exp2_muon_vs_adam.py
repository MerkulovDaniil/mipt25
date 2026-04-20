"""
Experiment 2: Muon vs Adam vs SGD+M on small MLP training
2-layer MLP on synthetic classification, numpy implementation.
Shows that Muon converges faster due to orthogonalized matrix updates.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_moons

np.random.seed(42)

# Dataset: 2-moons
X_train, y_train = make_moons(n_samples=2000, noise=0.15, random_state=42)
X_test, y_test = make_moons(n_samples=500, noise=0.15, random_state=123)

# Normalize
X_mean, X_std = X_train.mean(0), X_train.std(0)
X_train = (X_train - X_mean) / X_std
X_test = (X_test - X_mean) / X_std

# One-hot encode
n_classes = 2
Y_train = np.eye(n_classes)[y_train]
Y_test = np.eye(n_classes)[y_test]

n_train = len(X_train)
input_dim = 2
hidden_dim = 64
output_dim = 2

def relu(x):
    return np.maximum(0, x)

def relu_grad(x):
    return (x > 0).astype(float)

def softmax(x):
    e = np.exp(x - x.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)

def cross_entropy(probs, labels):
    return -np.mean(np.sum(labels * np.log(probs + 1e-12), axis=1))

def init_weights(seed=42):
    rng = np.random.RandomState(seed)
    # Xavier initialization
    W1 = rng.randn(input_dim, hidden_dim) * np.sqrt(2.0 / input_dim)
    b1 = np.zeros(hidden_dim)
    W2 = rng.randn(hidden_dim, output_dim) * np.sqrt(2.0 / hidden_dim)
    b2 = np.zeros(output_dim)
    return [W1, b1, W2, b2]

def forward(params, X):
    W1, b1, W2, b2 = params
    z1 = X @ W1 + b1
    a1 = relu(z1)
    z2 = a1 @ W2 + b2
    probs = softmax(z2)
    return probs, (z1, a1, z2)

def backward(params, X, Y, cache):
    W1, b1, W2, b2 = params
    z1, a1, z2 = cache
    batch = X.shape[0]

    probs = softmax(z2)
    dz2 = (probs - Y) / batch
    dW2 = a1.T @ dz2
    db2 = dz2.sum(0)

    da1 = dz2 @ W2.T
    dz1 = da1 * relu_grad(z1)
    dW1 = X.T @ dz1
    db1 = dz1.sum(0)

    return [dW1, db1, dW2, db2]

def evaluate(params, X, Y):
    probs, _ = forward(params, X)
    loss = cross_entropy(probs, Y)
    acc = np.mean(np.argmax(probs, axis=1) == np.argmax(Y, axis=1))
    return loss, acc

# Newton-Schulz for Muon
def newton_schulz5(G, num_iters=5):
    a, b, c = (3.4445, -4.7750, 2.0315)
    X = G.copy()
    nrm = np.linalg.norm(X, 'fro')
    if nrm < 1e-12:
        return X
    X = X / nrm
    X = X / (max(X.shape)**0.5)
    for _ in range(num_iters):
        A_ = X @ X.T
        B_ = b * A_ + c * A_ @ A_
        X = a * X + B_ @ X
    return X

n_steps = 4000
batch_size = 64
record_every = 20
n_seeds = 3

def train_adam(lr, seed=0):
    params = init_weights(seed)
    ms = [np.zeros_like(p) for p in params]
    vs = [np.zeros_like(p) for p in params]
    beta1, beta2, eps = 0.9, 0.999, 1e-8
    rng = np.random.RandomState(seed + 100)

    train_losses, test_losses = [], []
    for t in range(n_steps):
        if t % record_every == 0:
            tl, _ = evaluate(params, X_train, Y_train)
            tel, _ = evaluate(params, X_test, Y_test)
            train_losses.append(tl)
            test_losses.append(tel)

        idx = rng.choice(n_train, batch_size, replace=False)
        _, cache = forward(params, X_train[idx])
        grads = backward(params, X_train[idx], Y_train[idx], cache)

        for i in range(len(params)):
            ms[i] = beta1 * ms[i] + (1 - beta1) * grads[i]
            vs[i] = beta2 * vs[i] + (1 - beta2) * grads[i]**2
            mh = ms[i] / (1 - beta1**(t+1))
            vh = vs[i] / (1 - beta2**(t+1))
            params[i] = params[i] - lr * mh / (np.sqrt(vh) + eps)

    return train_losses, test_losses

def train_sgdm(lr, mom=0.9, seed=0):
    params = init_weights(seed)
    velocities = [np.zeros_like(p) for p in params]
    rng = np.random.RandomState(seed + 100)

    train_losses, test_losses = [], []
    for t in range(n_steps):
        if t % record_every == 0:
            tl, _ = evaluate(params, X_train, Y_train)
            tel, _ = evaluate(params, X_test, Y_test)
            train_losses.append(tl)
            test_losses.append(tel)

        idx = rng.choice(n_train, batch_size, replace=False)
        _, cache = forward(params, X_train[idx])
        grads = backward(params, X_train[idx], Y_train[idx], cache)

        for i in range(len(params)):
            velocities[i] = mom * velocities[i] + grads[i]
            params[i] = params[i] - lr * velocities[i]

    return train_losses, test_losses

def train_muon(lr_matrix, lr_vector, mom=0.95, seed=0):
    """Muon for matrix params, Adam for vector params (biases)."""
    params = init_weights(seed)
    # Muon state for matrices (W1, W2 at indices 0, 2)
    bufs = [np.zeros_like(p) for p in params]
    # Adam state for biases (b1, b2 at indices 1, 3)
    ms = [np.zeros_like(p) for p in params]
    vs = [np.zeros_like(p) for p in params]
    beta1, beta2, eps = 0.9, 0.999, 1e-8
    rng = np.random.RandomState(seed + 100)

    matrix_indices = [0, 2]  # W1, W2
    vector_indices = [1, 3]  # b1, b2

    train_losses, test_losses = [], []
    for t in range(n_steps):
        if t % record_every == 0:
            tl, _ = evaluate(params, X_train, Y_train)
            tel, _ = evaluate(params, X_test, Y_test)
            train_losses.append(tl)
            test_losses.append(tel)

        idx = rng.choice(n_train, batch_size, replace=False)
        _, cache = forward(params, X_train[idx])
        grads = backward(params, X_train[idx], Y_train[idx], cache)

        for i in range(len(params)):
            if i in matrix_indices:
                # Muon update
                bufs[i] = mom * bufs[i] + grads[i]
                G_tilde = mom * bufs[i] + grads[i]
                O = newton_schulz5(G_tilde)
                params[i] = params[i] - lr_matrix * O
            else:
                # Adam for biases
                ms[i] = beta1 * ms[i] + (1 - beta1) * grads[i]
                vs[i] = beta2 * vs[i] + (1 - beta2) * grads[i]**2
                mh = ms[i] / (1 - beta1**(t+1))
                vh = vs[i] / (1 - beta2**(t+1))
                params[i] = params[i] - lr_vector * mh / (np.sqrt(vh) + eps)

    return train_losses, test_losses

# Grid search for best lr
print("Tuning Adam...")
best_adam = (float('inf'), None)
for lr in [0.0003, 0.001, 0.003, 0.01]:
    tl, _ = train_adam(lr, seed=0)
    final = tl[-1]
    print(f"  Adam lr={lr}: {final:.4f}")
    if final < best_adam[0]:
        best_adam = (final, lr)

print("Tuning SGD+M...")
best_sgdm = (float('inf'), None)
for lr in [0.001, 0.003, 0.01, 0.03, 0.1]:
    tl, _ = train_sgdm(lr, seed=0)
    final = tl[-1]
    print(f"  SGD+M lr={lr}: {final:.4f}")
    if final < best_sgdm[0]:
        best_sgdm = (final, lr)

print("Tuning Muon...")
best_muon = (float('inf'), None, None)
for lr_m in [0.001, 0.003, 0.01, 0.03, 0.1]:
    lr_v = 0.001  # Adam lr for biases
    tl, _ = train_muon(lr_m, lr_v, seed=0)
    final = tl[-1]
    print(f"  Muon lr_m={lr_m}: {final:.4f}")
    if final < best_muon[0]:
        best_muon = (final, lr_m, lr_v)

print(f"\nBest: Adam lr={best_adam[1]}, SGD+M lr={best_sgdm[1]}, Muon lr_m={best_muon[1]}")

# Run with best lr, average over seeds
all_adam, all_sgdm, all_muon = [], [], []
all_adam_test, all_sgdm_test, all_muon_test = [], [], []
for seed in range(n_seeds):
    tl, tel = train_adam(best_adam[1], seed)
    all_adam.append(tl); all_adam_test.append(tel)
    tl, tel = train_sgdm(best_sgdm[1], seed=seed)
    all_sgdm.append(tl); all_sgdm_test.append(tel)
    tl, tel = train_muon(best_muon[1], best_muon[2], seed=seed)
    all_muon.append(tl); all_muon_test.append(tel)

adam_med = np.median(all_adam, axis=0)
sgdm_med = np.median(all_sgdm, axis=0)
muon_med = np.median(all_muon, axis=0)
adam_test = np.median(all_adam_test, axis=0)
sgdm_test = np.median(all_sgdm_test, axis=0)
muon_test = np.median(all_muon_test, axis=0)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

steps = np.arange(0, n_steps, record_every)

# Train loss
ax1.plot(steps, sgdm_med, '-', color='#e67e22', linewidth=2.5, label=f'SGD+Momentum (lr={best_sgdm[1]})', alpha=0.85)
ax1.plot(steps, adam_med, '-', color='#27ae60', linewidth=2.5, label=f'Adam (lr={best_adam[1]})', alpha=0.85)
ax1.plot(steps, muon_med, '-', color='#2980b9', linewidth=2.5, label=f'Muon (lr={best_muon[1]})', alpha=0.85)
ax1.set_xlabel('Steps', fontsize=14)
ax1.set_ylabel('Train Loss (CE)', fontsize=14)
ax1.set_title('Training Loss', fontsize=14)
ax1.legend(fontsize=11, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.tick_params(labelsize=12)

# Test loss
ax2.plot(steps, sgdm_test, '-', color='#e67e22', linewidth=2.5, label=f'SGD+Momentum', alpha=0.85)
ax2.plot(steps, adam_test, '-', color='#27ae60', linewidth=2.5, label=f'Adam', alpha=0.85)
ax2.plot(steps, muon_test, '-', color='#2980b9', linewidth=2.5, label=f'Muon', alpha=0.85)
ax2.set_xlabel('Steps', fontsize=14)
ax2.set_ylabel('Test Loss (CE)', fontsize=14)
ax2.set_title('Test Loss', fontsize=14)
ax2.legend(fontsize=11, loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.tick_params(labelsize=12)

fig.suptitle('2-layer MLP (hidden=64) on make_moons', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('/root/mipt25_work/files/exp_muon_vs_adam.pdf', bbox_inches='tight', dpi=150)
plt.savefig('/root/mipt25_work/files/exp_muon_vs_adam.png', bbox_inches='tight', dpi=150)
print("Done!")
