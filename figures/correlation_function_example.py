import numpy as np
import matplotlib.pyplot as plt

k = np.linspace(0, 300, 500)  # relative momentum in MeV/c

# Attractive: peak above 1 at low k, then approaches 1
def C_attractive(k):
    return 1 + 0.6 * np.exp(-k**2 / 800)

# Repulsive: dip below 1 at low k, then approaches 1
def C_repulsive(k):
    return 1 - 0.4 * np.exp(-k**2 / 1200)

fig, ax = plt.subplots(figsize=(4, 4))

ax.plot(k, C_attractive(k), label="Attractive pair", linewidth=2, color="steelblue")
ax.plot(k, C_repulsive(k), label="Repulsive pair", linewidth=3, color="tomato")
ax.axhline(1, color="gray", linestyle="--", linewidth=0.8, label="$\mathit{C}(\mathit{k}*)$ = 1")

ax.set_xlabel("$\mathit{k}$* (MeV/$\mathit{c}$)", fontsize=13)
ax.set_ylabel("$\mathit{C}(\mathit{k}*)$", fontsize=13)
ax.set_title("Femtoscopic Momentum Correlation Function", fontsize=14)
ax.set_yticks([])
ax.set_xlim(k[0], k[-1])
ax.legend()
plt.tight_layout()
plt.savefig("final_plots/femtoscopy_correlation.pdf", dpi=150)
plt.show()