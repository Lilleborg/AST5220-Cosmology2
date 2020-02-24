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
x_array_ticks = np.append(np.linspace(x_array.min(), x_array.max(), 5),0)
z_array = 1/np.exp(x_array) - 1
z_non_neg = z_array[z_array > 0]
a_array = np.exp(x_array)
km = 1e3                                # [m]
Mpc = 3.08567758e22                     # [m]
h = 0.7                                 # unitless hubble
H0 = 100 * h * km/Mpc                   # Hubble parameter today
c = 299792458                           # Speed of light [m/s]
Omega_matter = OmegaB + OmegaCDM        # Omega matter combined
# Some important indices and scaling of x_array between 0 and 1:
idx_rad_mat_eq = np.argmin(np.abs(OmegaR-Omega_matter))
idx_mat_lambda_eq = np.argmin(np.abs(OmegaLambda-Omega_matter))
delta_x = x_array.max() - x_array.min()
xmax_rad_scaled = (x_array[idx_rad_mat_eq]+np.abs(x_array.min()))/delta_x
xmax_mat_scaled = (x_array[idx_mat_lambda_eq]+np.abs(x_array.min()))/delta_x

# Benchmarking/testing
# --------------------

# H(a_0) == H0
idx_a0 = np.argmin(np.abs(a_array-1))   # closest index to a == 1
print('-'*29)
print('{:^29s}'.format('H(a0) == H0'))
print('-'*29)
print('Closest a = a0:   {:.5e}'.format(a_array[idx_a0]))
print('True value of H0: {:.5e}'.format(H0))
print('Closest H(a0):    {:.5e}'.format(H_of_x[idx_a0]))
print('Ratio H(a0)/H0:   {:.5e}'.format(H_of_x[idx_a0]/H0))

# a*H/Hp == 1
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(x_array, Hp_of_x/(a_array*H_of_x), linewidth=0.7)
ax.set_xlabel(r'$x = ln(a)$')
ax.set_ylabel(r'$\mathcal{H}/(aH)$')
ax.tick_params(axis='y', labelrotation=45)
ax.set_xticks(x_array_ticks)

fig.savefig('./figs/ratio_Hprime_aH.pdf')

# Ratio of Hp and its derivatives in each regime
fig, axes = plt.subplots(2, figsize=(6, 9))
fig.suptitle(r'Ratio of $\mathcal{H}$ and its derivatives')

axes[0].set_title(r'$(1-N/2)$ in each regime')
axes[0].plot(x_array, dHpdx_of_x/Hp_of_x, color='C2')
axes[0].set_ylabel(r'$\mathcal{H}^{\prime}/\mathcal{H}$')
axes[0].axhline(y=-1, label=r'$N=4$', color='C0', ls='--', alpha=0.7, xmax=xmax_rad_scaled)
axes[0].axhline(y=-0.5, label=r'$N=3$', color='C1', ls='--', alpha=0.7, xmin=xmax_rad_scaled, xmax=xmax_mat_scaled)
axes[0].axhline(y=1, label=r'$N=0$', color='C3', ls='--', alpha=0.7, xmin=xmax_mat_scaled)

axes[1].set_title(r'$(1-N/2)^2$ in each regime')
axes[1].plot(x_array, ddHpddx_of_x/Hp_of_x, color='C2')
axes[1].set_ylabel(r'$\mathcal{H}^{\prime\prime}/\mathcal{H}$')
axes[1].axhline(y=1, label=r'$N=4$', color='C0', ls='--', alpha=0.7, xmax=xmax_rad_scaled)
axes[1].axhline(y=0.25, label=r'$N=3$', color='C1', ls='--', alpha=0.7, xmin=xmax_rad_scaled, xmax=xmax_mat_scaled)
axes[1].axhline(y=1, label=r'$N=0$', color='C3', ls='-.', alpha=0.7, xmin=xmax_mat_scaled)

for ax in axes.flatten():
    ax.axvspan(x_array.min(), x_array[idx_rad_mat_eq], alpha=0.25, color='C0')
    ax.axvspan(x_array[idx_rad_mat_eq], x_array[idx_mat_lambda_eq], alpha=0.25, color='C1')
    ax.axvspan(x_array[idx_mat_lambda_eq], x_array.max(), alpha=0.25, color='C3')
    ax.set_xticks(x_array_ticks)
    ax.margins(x=0)
    ax.legend()

fig.savefig('./figs/ratio_Hprime_and_derivatives.pdf')

# Ratio eta*Hp/c = 1 back in time
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(x_array, eta_of_x*Hp_of_x/c, color='C2')
ax.axhline(y=1, color='C0', ls='--', alpha=0.7, label='Convergence value')
ax.set_ylabel(r'$\eta \mathcal{H} /c$')
ax.set_xlabel(r'$x = ln(a)$')
ax.set_yscale('log')
ax.legend()

fig.savefig('./figs/ratio_eta_Hp_over_c.pdf')

plt.show()
