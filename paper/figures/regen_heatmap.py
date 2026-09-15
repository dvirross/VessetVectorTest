"""Regenerate fig8_chisq_heatmap with § instead of † for analytically predicted pairs."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PATTERN_LABELS = ['Haflaga', 'Dilug', 'Week', 'Week-\nDilug', 'DiD']
PATTERN_LABELS_CLEAN = ['Haflaga', 'Dilug', 'Week', 'Week-Dilug', 'DiD']

# Spearman correlation data (rate-based, seed 17, B=50,000)
# From paper Table tab:chisq
pairs_data = [
    ('Haflaga', 'Dilug-in-Dilug', 0.27, 0.003, 0.030, True),
    ('Haflaga', 'Week-Dilug',     0.23, 0.010, 0.090, True),
    ('Dilug',   'Week-Dilug',    -0.12, 0.206, 1.000, False),
    ('Haflaga', 'Week',           0.10, 0.272, 1.000, False),
    ('Dilug',   'Week',           0.09, 0.325, 1.000, False),
    ('Week',    'Dilug-in-Dilug', 0.09, 0.326, 1.000, False),
    ('Week',    'Week-Dilug',     0.07, 0.474, 1.000, False),
    ('Week-Dilug', 'Dilug-in-Dilug', 0.07, 0.483, 1.000, False),
    ('Haflaga', 'Dilug',          0.04, 0.676, 1.000, False),
    ('Dilug',   'Dilug-in-Dilug',-0.02, 0.815, 1.000, False),
]

n = 5
idx = {
    'Haflaga': 0, 'Dilug': 1, 'Week': 2,
    'Week-Dilug': 3, 'Dilug-in-Dilug': 4
}

rs_mat   = np.zeros((n, n))
holm_mat = np.ones((n, n))
pred_mat = np.zeros((n, n), dtype=bool)

for a, b, rs, p_perm, p_holm, predicted in pairs_data:
    i, j = idx[a], idx[b]
    rs_mat[i, j] = rs_mat[j, i] = rs
    holm_mat[i, j] = holm_mat[j, i] = p_holm
    if predicted:
        pred_mat[i, j] = pred_mat[j, i] = True

# Axis labels (split Week-Dilug for display)
ax_labels = ['Haflaga', 'Dilug', 'Week', 'Week-\nDilug', 'DiD']

fig, ax = plt.subplots(figsize=(7, 6))
cmap = plt.cm.RdBu_r
im = ax.imshow(rs_mat, cmap=cmap, vmin=-0.35, vmax=0.35, aspect='auto')

ax.set_xticks(range(n))
ax.set_xticklabels(ax_labels, fontsize=10)
ax.set_yticks(range(n))
ax.set_yticklabels(ax_labels, fontsize=10)

# Draw grid lines
for k in range(n + 1):
    ax.axhline(k - 0.5, color='white', lw=1.5)
    ax.axvline(k - 0.5, color='white', lw=1.5)

# Annotate cells
for i in range(n):
    for j in range(n):
        if i == j:
            ax.text(j, i, '—', ha='center', va='center', fontsize=11,
                    color='#555555', fontweight='bold')
            continue
        rs = rs_mat[i, j]
        p_holm = holm_mat[i, j]
        is_pred = pred_mat[i, j]

        # Build label: sign + value + § (if predicted) + * (if Holm sig)
        sign = '+' if rs >= 0 else ''
        label = f'{sign}{rs:.2f}'
        if is_pred:
            label += '§'   # § section sign
        if p_holm < 0.05:
            label += '*'

        # Text colour: white on dark cells, black on light
        bg_intensity = abs(rs) / 0.35
        color = 'white' if bg_intensity > 0.55 else 'black'
        weight = 'bold' if p_holm < 0.05 else 'normal'
        ax.text(j, i, label, ha='center', va='center',
                fontsize=9, color=color, fontweight=weight)

cb = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cb.set_label(r'$r_s$', fontsize=11, rotation=0, labelpad=10)

ax.set_title('Pairwise Spearman correlations (rate-based)\n'
             '§ analytically predicted * Holm adj $p < .05$',
             fontsize=11, pad=10)

plt.tight_layout()
plt.savefig('/home/user/VessetVectorTest/paper/figures/fig8_chisq_heatmap.pdf',
            bbox_inches='tight')
plt.savefig('/home/user/VessetVectorTest/paper/figures/fig8_chisq_heatmap.png',
            dpi=150, bbox_inches='tight')
print('Saved.')
