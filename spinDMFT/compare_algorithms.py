import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Mapping from beta to (reference_file, mh_file)
files = {
    1: ("Data/spinmodel=ISO__beta=1.hdf5", "Data/spinmodel=ISO__beta=1XXX.hdf5"),
    3: ("Data/spinmodel=ISO__beta=3.hdf5", "Data/spinmodel=ISO__beta=3XX.hdf5"),
    5: (None, "Data/spinmodel=ISO__beta=5XXX.hdf5"),  # No reference at beta=5
}

fig, axes = plt.subplots(2, 3, figsize=(13, 7))
fig.suptitle('Correlation Comparison: MH vs Original (zz-component)', fontsize=14)

for idx, (beta, (ref_file, mh_file)) in enumerate(files.items()):
    ax = axes.flat[idx]
    
    # Read MH data
    if Path(mh_file).exists():
        with h5py.File(mh_file, 'r') as h:
            corr_mh = h['results/Re_correlation'][()]
            num_ts = h['parameters/num_TimePoints'][()]
            tau = np.arange(num_ts) * h['parameters/delta_t'][()]
            # Extract zz component (index 2)
            corr_mh_zz = corr_mh[2::3] if len(corr_mh.shape)==1 else corr_mh[2,:]
        
        ax.plot(tau/beta, corr_mh_zz, 'o-', label='MH', color='red', markersize=4, linewidth=1.5, alpha=0.7)
        ax.set_xlabel(r'$\tau/\beta$')
        ax.set_ylabel(r'$g^{zz}(\tau)$')
        ax.set_title(f'β={beta}', fontsize=11)
        ax.grid(True, alpha=0.3)
    
    # Read reference data if available
    if ref_file and Path(ref_file).exists():
        with h5py.File(ref_file, 'r') as h:
            corr_ref = h['results/Re_correlation'][()]
            num_ts_ref = h['parameters/num_TimePoints'][()]
            tau_ref = np.arange(num_ts_ref) * h['parameters/delta_t'][()]
            corr_ref_zz = corr_ref[2::3] if len(corr_ref.shape)==1 else corr_ref[2,:]
        
        ax.plot(tau_ref/beta, corr_ref_zz, 's--', label='Original (IS)', color='blue', markersize=4, linewidth=1.5, alpha=0.7)
        ax.legend(fontsize=9)
    elif not ref_file:
        ax.text(0.5, 0.5, 'No reference data', transform=ax.transAxes, ha='center', va='center', fontsize=10, color='gray')

plt.tight_layout()
plt.savefig('Plots/mh_vs_original.pdf', dpi=300)
print("✓ Saved comparison plot to Plots/mh_vs_original.pdf")

# Print summary stats
print("\nSummary:")
for beta, (ref_file, mh_file) in files.items():
    if Path(mh_file).exists():
        with h5py.File(mh_file, 'r') as h:
            print(f"\nbeta={beta} (MH):")
            print(f"  num_Samples={h['parameters/num_Samples'][()]}")
            print(f"  num_TimeSteps={h['parameters/num_TimeSteps'][()]}")
            try:
                acc = h['results/acceptance_rate'][()]
                print(f"  acceptance_rate={acc:.3f}")
            except:
                pass

