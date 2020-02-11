import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Plot styling
plt.style.use('bmh')
mpl.rc('font', family='serif', size=15)
mpl.rc('text', usetex=False)

# Load data
# ---------
x_array, eta_of_x, H_of_x, dHdx_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, OmegaB, OmegaCDM, \
    OmegaLambda, OmegaR, _, _ = np.loadtxt("../data/cosmology.txt", unpack=True)

# Data handling and converting
# ----------------------------
x_array_ticks = np.linspace(x_array.min(), x_array.max(), 6)
z_array = 1/np.exp(x_array) - 1
z_non_neg = z_array[z_array > 0]
a_array = np.exp(x_array)
km = 1e3                                # [m]
Mpc = 3.08567758e22                     # [m]
# H-prime in presentable units [km/s/Mpc]:
Hp_units = Hp_of_x * Mpc/km
H_units = Hp_units / a_array            # H (not prime)  [km/s/Mpc]
eta_units = eta_of_x / (Mpc*1e3)        # conformal time [Gpc]
# Indices for radiation-matter and matter-lambda equalities:
Omega_matter = OmegaB + OmegaCDM        # Omega matter combined
idx_rad_mat_eq = np.argmin(np.abs(OmegaR-Omega_matter))
idx_mat_lambda_eq = np.argmin(np.abs(OmegaLambda-Omega_matter))

# Plotting
# --------
# Omegas vs x
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(x_array, OmegaR, label=r'$\Omega_R$')
ax.plot(x_array, OmegaCDM, label=r'$\Omega_{CDM}$')
ax.plot(x_array, OmegaB, label=r'$\Omega_B$')
ax.plot(x_array, OmegaLambda, label=r'$\Omega_\Lambda$')
ax.plot(x_array, Omega_matter, label=r'$\Omega_{\rm{Matter}}$', ls='--', alpha=0.7)
ax.axvspan(x_array.min(), x_array[idx_rad_mat_eq], alpha=0.25, color='C0')
ax.axvspan(x_array[idx_rad_mat_eq], x_array[idx_mat_lambda_eq], alpha=0.25, color='C1')
ax.axvspan(x_array[idx_mat_lambda_eq], x_array.max(), alpha=0.25, color='C3')
ax.legend(loc=6)
ax.set_xlabel(r'$x = ln(a)$')
ax.set_ylabel(r'$\Omega_i = \frac{\rho_i}{\rho_c}$')
ax.set_xticks(x_array_ticks)

fig.savefig('./figs/omegas_of_x.pdf', bbox_inces='tight')

# Plotting Hubble parameter vs x
fig, axes = plt.subplots(2, 2)
axes[0, 0].plot(x_array, H_units, color='C2')
axes[0, 0].set_ylabel(r'$H(x) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')
axes[0, 0].set_xticks(x_array_ticks)

# Plotting Hubble parameter vs z
axes[0, 1].plot(z_non_neg, H_units[z_array > 0], color='C2')
axes[0, 1].set_xlabel(r'$z = a^{-1} - 1$')
axes[0, 1].set_ylabel(r'$H(z) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')
axes[0, 1].set_xlim(z_non_neg.max()+z_non_neg.max()*3, z_non_neg.min()+0.0001)
axes[0, 1].set_xscale('log')

# Plotting conformal time
axes[1, 0].plot(x_array, eta_units, color='C2')
axes[1, 0].set_ylabel(r'$\eta(x) [\rm{Gpc}]$')
axes[1, 0].set_xticks(x_array_ticks)

# PLotting H-prime
axes[1, 1].plot(x_array, Hp_units, color='C2')
axes[1, 1].set_ylabel(r'$\mathcal{H}(x) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')
axes[1, 1].set_xticks(x_array_ticks)

# Filling background in each axis
for i, ax in enumerate(axes.flatten()):
    ax.set_yscale('log')
    if i != 1:
        ax.axvspan(x_array.min(), x_array[idx_rad_mat_eq], alpha=0.25, color='C0')
        ax.axvspan(x_array[idx_rad_mat_eq], x_array[idx_mat_lambda_eq], alpha=0.25, color='C1')
        ax.axvspan(x_array[idx_mat_lambda_eq], x_array.max(), alpha=0.25, color='C3')
        ax.set_xlabel(r'$x = ln(a)$')
    else:
        ax.axvspan(z_array.max(), z_array[idx_rad_mat_eq], alpha=0.25, color='C0')
        ax.axvspan(z_array[idx_rad_mat_eq], z_array[idx_mat_lambda_eq], alpha=0.25, color='C1')
        ax.axvspan(z_array[idx_mat_lambda_eq], z_array.min(), alpha=0.25, color='C3')

fig.savefig('./figs/Hubble_eta_of_x.pdf', bbox_inces='tight')
plt.show()
