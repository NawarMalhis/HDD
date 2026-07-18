#!/usr/bin/env python3
"""
Binding site length distribution in training datasets
- Very short (<=30 AA)
- Medium short (31-70 AA)
- Long (>70 AA)
"""

import matplotlib.pyplot as plt
import numpy as np

# ========================== DATA (example / placeholder) ==========================
# Replace with your actual counts from the annotated datasets (e.g. via aff_load3)
datasets = ["DisProt", "CAID1u", "CAID23u", "DBs", "TR2008u", "TR2017"]

# Percentages of binding residues in each category
very_short = [7, 12, 15, 10, 100, 100]  # [100, 100, 10, 15, 12]   # <= 30 AA
medium_short = [12, 9, 8, 10, 0, 0]  # [0, 0, 10, 8, 9]       # 31–70 AA
long_sites = [81, 79, 77, 80, 0, 0]  # [0, 0, 80, 77, 79]       # >70 AA

# Verify sums ≈ 100%
for i in range(len(datasets)):
    total = very_short[i] + medium_short[i] + long_sites[i]
    print(f"{datasets[i]:8} → {total:3.0f}%")

# ========================== PLOTTING ==========================
fig, ax = plt.subplots(figsize=(9, 6))

x = np.arange(len(datasets))
width = 0.75

# Stacked bars
p1 = ax.bar(x, very_short, width, label='Very short (≤30 AA)', color='#1f77b4')
p2 = ax.bar(x, medium_short, width, bottom=very_short,
            label='Medium short (31–70 AA)', color='#ff7f0e')
p3 = ax.bar(x, long_sites, width,
            bottom=np.array(very_short) + np.array(medium_short),
            label='Long (>70 AA)', color='#2ca02c')

ax.set_ylabel("% of binding residues", fontsize=14)
ax.set_title("(training datasets)", fontsize=15, pad=20)
ax.set_xticks(x)
ax.set_xticklabels(datasets, fontsize=12, rotation=0)
ax.set_ylim(0, 105)
ax.legend(loc='upper center', fontsize=11)
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Add percentage labels on bars
for i in range(len(datasets)):
    cum = 0
    for perc, color in zip([very_short[i], medium_short[i], long_sites[i]],
                           ['#1f77b4', '#ff7f0e', '#2ca02c']):
        if perc > 3:  # only label larger segments
            ax.text(i, cum + perc/2, f"{perc}%", ha='center', va='center',
                    fontsize=10, color='white', fontweight='bold')
        cum += perc

plt.tight_layout()
plt.savefig("Figure_3_left_training_site_lengths.png", dpi=300, bbox_inches='tight')
plt.show()