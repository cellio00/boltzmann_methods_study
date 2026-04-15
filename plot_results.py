#!/usr/bin/env python3
"""
Enhanced visualization for two-beam neutrino collision
Shows right-moving electron beam and left-moving x-flavor beam
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

print("Loading output data...")
data = np.loadtxt("output.dat")

# Parse columns
it = data[:, 0].astype(int)
t = data[:, 1]
ix = data[:, 2].astype(int)
x = data[:, 3]

# Distribution function components
f_ee_plus = data[:, 4]    # Electron, moving right
f_ex_plus = data[:, 5]
f_ex_conj_plus = data[:, 6]
f_xx_plus = data[:, 7]    # X-flavor, moving right

f_ee_minus = data[:, 8]   # Electron, moving left
f_ex_minus = data[:, 9]
f_ex_conj_minus = data[:, 10]
f_xx_minus = data[:, 11]  # X-flavor, moving left

# Moments
J_e = data[:, 12]
J_x = data[:, 15]

# Get timesteps
timesteps = np.unique(it)
n_times = len(timesteps)

print(f"Data loaded: {n_times} timesteps")

# ==============================================================================
# Plot 1: Two-Beam Visualization at Final Time
# ==============================================================================

final_it = timesteps[-1]
mask = (it == final_it)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f'Two-Beam Neutrino Collision (t = {t[mask][0]:.3f})', fontsize=16, fontweight='bold')

# Panel 1: Right-moving beam (f_plus)
ax = axes[0, 0]
# Original electron beam (started on left)
ax.plot(x[mask], f_ee_plus[mask], 'b-', linewidth=2.5, label='$f_{ee}^+$ (original e beam →)')
# X-flavor created by oscillations
ax.plot(x[mask], f_xx_plus[mask], 'r-', linewidth=2.5, alpha=0.7, label='$f_{xx}^+$ (oscillated to x →)')
# Total right-moving
total_plus = f_ee_plus[mask] + f_xx_plus[mask]
ax.plot(x[mask], total_plus, 'k:', linewidth=2, label='Total (→)')
ax.set_xlabel('Position x', fontsize=12)
ax.set_ylabel('Distribution (right-moving)', fontsize=12)
ax.set_title('Right-Moving Neutrinos (f⁺)', fontsize=13, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.2)

# Panel 2: Left-moving beam (f_minus)
ax = axes[0, 1]
# X-flavor created by oscillations  
ax.plot(x[mask], f_ee_minus[mask], 'b-', linewidth=2.5, alpha=0.7, label='$f_{ee}^-$ (oscillated to e ←)')
# Original x-flavor beam (started on right)
ax.plot(x[mask], f_xx_minus[mask], 'r-', linewidth=2.5, label='$f_{xx}^-$ (original x beam ←)')
# Total left-moving
total_minus = f_ee_minus[mask] + f_xx_minus[mask]
ax.plot(x[mask], total_minus, 'k:', linewidth=2, label='Total (←)')
ax.set_xlabel('Position x', fontsize=12)
ax.set_ylabel('Distribution (left-moving)', fontsize=12)
ax.set_title('Left-Moving Neutrinos (f⁻)', fontsize=13, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.2)



# Panel 3: Flavor densities (moments)
ax = axes[1, 0]
ax.plot(x[mask], J_e[mask], 'b-', linewidth=3, label='$J_e$ (electron density)')
ax.plot(x[mask], J_x[mask], 'r-', linewidth=3, label='$J_x$ (x-flavor density)')
ax.plot(x[mask], J_e[mask] + J_x[mask], 'k:', linewidth=2, label='$J_{total}$ (conserved)')
ax.set_xlabel('Position x', fontsize=12)
ax.set_ylabel('Number Density', fontsize=12)
ax.set_title('Total Flavor Densities (both beams)', fontsize=13, fontweight='bold')
ax.legend(loc='upper right', fontsize=11)
ax.grid(True, alpha=0.3)
ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5, label='Initial boundary')

# Panel 4: Flavor mixing (oscillation signature)
ax = axes[1, 1]

# Compute flavor purity: how much electron vs x-flavor
total = J_e[mask] + J_x[mask]
electron_fraction = np.divide(J_e[mask], total, where=total>1e-10, out=np.zeros_like(J_e[mask]))

ax.fill_between(x[mask], 0, electron_fraction, alpha=0.3, color='blue', label='Electron fraction')
ax.fill_between(x[mask], electron_fraction, 1, alpha=0.3, color='red', label='X-flavor fraction')
ax.plot(x[mask], electron_fraction, 'b-', linewidth=2)
ax.set_xlabel('Position x', fontsize=12)
ax.set_ylabel('Electron Fraction', fontsize=12)
ax.set_title('Flavor Composition', fontsize=13, fontweight='bold')
ax.set_ylim(-0.05, 1.05)
ax.legend(loc='right', fontsize=11)
ax.grid(True, alpha=0.3)
ax.axhline(y=0.5, color='purple', linestyle=':', alpha=0.7, linewidth=1.5, label='Equal mixing')

plt.tight_layout()
plt.savefig('two_beam_collision.png', dpi=150, bbox_inches='tight')
print("Saved: two_beam_collision.png")

# ==============================================================================
# Plot 2: Time Evolution Sequence
# ==============================================================================

# Select key timepoints to show evolution
n_snapshots = min(6, n_times)
snapshot_indices = np.linspace(0, n_times-1, n_snapshots, dtype=int)

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle('Two-Beam Evolution Sequence', fontsize=16, fontweight='bold')

for idx, snap_idx in enumerate(snapshot_indices):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]
    
    step = timesteps[snap_idx]
    mask = (it == step)
    
    # Plot the primary beams and their oscillation products
    # Right-moving: electron beam + oscillated x-flavor
    ax.plot(x[mask], f_ee_plus[mask], 'b-', linewidth=2.5, label='$f_{ee}^+$ (e→)')
    ax.plot(x[mask], f_xx_plus[mask], 'r--', linewidth=1.5, alpha=0.6, label='$f_{xx}^+$ (x→)')
    # Left-moving: x-flavor beam + oscillated electron
    ax.plot(x[mask], f_xx_minus[mask], 'r-', linewidth=2.5, label='$f_{xx}^-$ (x←)')
    ax.plot(x[mask], f_ee_minus[mask], 'b--', linewidth=1.5, alpha=0.6, label='$f_{ee}^-$ (e←)')
    
    ax.set_ylim(-0.05, 1.2)
    ax.set_xlabel('Position x', fontsize=9)
    ax.set_ylabel('Distribution', fontsize=9)
    ax.set_title(f't = {t[mask][0]:.3f}', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.3)
    
    if idx == 0:
        ax.legend(fontsize=7, loc='upper left', ncol=2)

plt.tight_layout()
plt.savefig('beam_evolution_sequence.png', dpi=150, bbox_inches='tight')
print("Saved: beam_evolution_sequence.png")

# ==============================================================================
# Plot 3: Animation
# ==============================================================================

print("\nCreating animation...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

def update(frame):
    step = timesteps[frame]
    mask = (it == step)
    
    for ax in axes:
        ax.clear()
    
    # Left panel: The two beams with oscillation products
    ax = axes[0]
    # Right-moving beam (started as electron on left)
    ax.plot(x[mask], f_ee_plus[mask], 'b-', linewidth=3, label='Electron beam → (original)')
    ax.plot(x[mask], f_xx_plus[mask], 'r--', linewidth=2, alpha=0.7, label='X-flavor → (from oscillation)')
    # Left-moving beam (started as x-flavor on right)
    ax.plot(x[mask], f_xx_minus[mask], 'r-', linewidth=3, label='X-flavor beam ← (original)')
    ax.plot(x[mask], f_ee_minus[mask], 'b--', linewidth=2, alpha=0.7, label='Electron ← (from oscillation)')
    
    ax.set_ylim(-0.05, 1.3)
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Distribution', fontsize=12)
    ax.set_title(f'Two-Beam Collision (t = {t[mask][0]:.3f})', fontsize=13, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.3, linewidth=1)
    
    # Right panel: Total densities
    ax = axes[1]
    ax.fill_between(x[mask], 0, J_e[mask], alpha=0.4, color='blue', label='Electron')
    ax.fill_between(x[mask], J_e[mask], J_e[mask]+J_x[mask], alpha=0.4, color='red', label='X-flavor')
    ax.plot(x[mask], J_e[mask] + J_x[mask], 'k-', linewidth=2, label='Total')
    ax.set_ylim(-0.05, 1.3)
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Number Density', fontsize=12)
    ax.set_title('Flavor Densities', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    return axes

anim = FuncAnimation(fig, update, frames=len(timesteps), interval=50, blit=False)

plt.tight_layout()

# Save animation
if len(sys.argv) > 1 and sys.argv[1] == '--save-gif':
    print("Saving animation as two_beam_collision.gif...")
    anim.save('two_beam_collision.gif', writer='pillow', fps=20)
    print("Saved: two_beam_collision.gif")
else:
    print("\nShowing animation (close window to exit)...")
    print("To save as GIF, run: python plot_two_beams.py --save-gif")
    plt.show()

print("\nDone!")
print("\nGenerated files:")
print("  - two_beam_collision.png       (final state)")
print("  - beam_evolution_sequence.png  (6 time snapshots)")
if len(sys.argv) > 1 and sys.argv[1] == '--save-gif':
    print("  - two_beam_collision.gif       (animation)")