import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
fig.suptitle('Correlation Comparison: MH vs Reference', fontsize=13, fontweight='bold')

files_to_compare = [
    (1, "Data/spinmodel=ISO__beta=1.hdf5", "Data/spinmodel=ISO__beta=1XXX.hdf5"),
    (3, "Data/spinmodel=ISO__beta=3.hdf5", "Data/spinmodel=ISO__beta=3XX.hdf5"),
    (5, None, "Data/spinmodel=ISO__beta=5XXX.hdf5"),
]

for plot_idx, (beta, ref_file, mh_file) in enumerate(files_to_compare):
    ax = axes[plot_idx]
    ax.set_title(f'β = {beta}', fontsize=11, fontweight='bold')
    ax.set_xlabel(r'Time index')
    ax.set_ylabel(r'Correlation (zz)')
    ax.grid(True, alpha=0.2)
    
    # Read and plot MH
    if Path(mh_file).exists():
        with h5py.File(mh_file, 'r') as h:
            corr = h['results/Re_correlation'][:]
            # Shape is (num_directions, num_TimePoints). Extract zz which is typically index 2 or 8
            # For a single spin: directions are xx, xy, xz, yx, yy, yz, zx, zy, zz
            # So zz should be index 8 if all 9 components are present, or index 0 if only zz
            if corr.shape[0] >= 9:
                corr_zz = corr[8, :]  # zz component
            else:
                corr_zz = corr[0, :]  # assume it's zz-only
            
            tau = np.arange(len(corr_zz))
            ax.plot(tau, corr_zz, 'o-', label='MH', color='red', markersize=6, linewidth=2, alpha=0.8)
    
    # Read and plot reference
    if ref_file and Path(ref_file).exists():
        with h5py.File(ref_file, 'r') as h:
            corr_ref = h['results/Re_correlation'][:]
            if corr_ref.shape[0] >= 9:
                corr_ref_zz = corr_ref[8, :]
            else:
                corr_ref_zz = corr_ref[0, :]
            
            tau_ref = np.arange(len(corr_ref_zz))
            ax.plot(tau_ref, corr_ref_zz, 's--', label='Reference (IS)', color='blue', markersize=5, linewidth=1.5, alpha=0.7)
    else:
        ax.text(0.5, 0.5, 'No reference', transform=ax.transAxes, ha='center', va='center', fontsize=10, color='gray')
    
    if ref_file:
        ax.legend(loc='best', fontsize=9)

plt.tight_layout()
plt.savefig('Plots/mh_comparison.pdf', dpi=300, bbox_inches='tight')
print("✓ Saved plot: Plots/mh_comparison.pdf")

# Print summary
print("\n" + "="*60)
print("MH BENCHMARK SUMMARY")
print("="*60)

for beta, ref_file, mh_file in files_to_compare:
    if Path(mh_file).exists():
        with h5py.File(mh_file, 'r') as h:
            corr = h['results/Re_correlation'][:]
            print(f"\nbeta = {beta}:")
            print(f"  Correlation shape: {corr.shape}")
            print(f"  Correlation values (zz, first 5 time points):")
            if corr.shape[0] >= 9:
                print(f"    {corr[8, :5]}")
            else:
                print(f"    {corr[0, :5]}")
    
    if ref_file and Path(ref_file).exists():
        with h5py.File(ref_file, 'r') as h:
            corr_ref = h['results/Re_correlation'][:]
            print(f"  Reference shape: {corr_ref.shape}")
            if corr_ref.shape[0] >= 9:
                print(f"  Reference (zz, first 5): {corr_ref[8, :5]}")
            else:
                print(f"  Reference (first 5): {corr_ref[0, :5]}")
