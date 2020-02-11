import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

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
h = 0.7                                 # unitless hubble
H0 = 100 * h * km/Mpc                   # Hubble parameter today
Omega_matter = OmegaB + OmegaCDM        # Omega matter combined
idx_rad_mat_eq = np.argmin(np.abs(OmegaR-Omega_matter))
idx_mat_lambda_eq = np.argmin(np.abs(OmegaLambda-Omega_matter))

# Benchmarking/testing
# --------------------

# H(a_0) == H0
idx_a0 = np.argmin(np.abs(a_array-1))   # closest index of a == 1
print('-'*29)
print('{:^29s}'.format('H(a0) == H0'))
print('-'*29)
print('Closest a = a0:   {:.5e}'.format(a_array[idx_a0]))
print('True value of H0: {:.5e}'.format(H0))
print('Closest H(a0):    {:.5e}'.format(H_of_x[idx_a0]))
print('Ratio H(a0)/H0:   {:.5e}'.format(H_of_x[idx_a0]/H0))

# a*H/Hp == 1
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(x_array, Hp_of_x/(a_array*H_of_x))
ax.set_xlabel(r'$x = ln(a)$')
ax.set_ylabel(r'$\mathcal{H}/(aH)$')
ax.tick_params(axis='y',labelrotation=45)

# Ratio of Hp and its derivatives in each regime
fig, axes = plt.subplots(2,figsize=(6, 9))
axes[0].plot(x_array, dHpdx_of_x/Hp_of_x, color='C2')
axes[0].set_ylabel(r'$\mathcal{H}^{\prime}/\mathcal{H}$')
axes[0].set_title(r'$(1-N)$ in each regime')
fig.suptitle(r'Ratio of $\mathcal{H}$ and its derivatives')

for ax in axes.flatten():
    ax.axhline(y=-1, label=r'$N=-2$', color='C0', ls='--')
    ax.axhline(y=-0.5, label=r'$N=-3/2$', color='C1', ls='--')
    ax.axhline(y=1, label=r'$N=0$', color='C3', ls='--')

    ax.axvspan(x_array.min(), x_array[idx_rad_mat_eq], alpha=0.25, color='C0')
    ax.axvspan(x_array[idx_rad_mat_eq], x_array[idx_mat_lambda_eq], alpha=0.25, color='C1')
    ax.axvspan(x_array[idx_mat_lambda_eq], x_array.max(), alpha=0.25, color='C3')

axes[0].legend()
plt.show()
