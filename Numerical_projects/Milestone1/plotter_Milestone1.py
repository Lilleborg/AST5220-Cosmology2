import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Plot styling
plt.style.use('bmh')
plt.rc('figure', figsize=(8,5.5))

# Load data
x_array, eta_of_x, Hp_of_x, _, OmegaB, OmegaCDM, OmegaLambda, OmegaR, _, _ = np.loadtxt(
    "../data/cosmology.txt", unpack=True)

# def redshift_to_x(z):
#     return np.log(1/(1+z))

# def x_to_redshift(x):
#     return 1/np.exp(x) - 1
#ax2 = ax.secondary_xaxis('top', functions=(x_to_redshift,redshift_to_x))

# Data handling and preparation for plotting
z_array = 1/np.exp(x_array) - 1
a_array = np.exp(x_array)

# Plotting omegas
fig, ax = plt.subplots()
ax.plot(x_array, OmegaB, label=r'$\Omega_b$')
ax.plot(x_array, OmegaCDM, label=r'$\Omega_{CDM}$')
ax.plot(x_array, OmegaLambda, label=r'$\Omega_\Lambda$')
ax.plot(x_array, OmegaR, label=r'$\Omega_R$')
ax.legend(loc=6)
ax.set_xlabel(r'x = ln(a)')

# Plotting Hubble parameter and conformal time
fig, ax = plt.subplots(2,2,figsize=(14,10))
ax[0,0].plot(x_array,Hp_of_x/a_array)
ax[0,0].set_xlabel(r'$x = ln(a)$')
ax[0,0].set_ylabel(r'$H(x)$')

ax[0,1].plot(z_array,Hp_of_x/a_array)
ax[0,1].set_xlim(z_array[-1],z_array[0])
ax[0,1].set_xlabel(r'$z = \frac{1}{a} - 1$')
ax[0,1].set_ylabel('$H(z)$')
ax[0,1].set_xlim(z_array.max()+z_array.max()*0.05,z_array.min()+0.0001-z_array.min()*0.05)
ax[0,1].set_xscale('log')


plt.show()


# ax2 = ax.twiny()
# ax2.plot(z_array, np.zeros(z_array.shape))
# ax2.plot(z_array, np.zeros(z_array.shape))
# ax2.plot(z_array, np.zeros(z_array.shape))
# ax2.plot(z_array, np.zeros(z_array.shape))
# ax2.set_xscale('log')
# #ax2.set_xlim(z_array[0],z_array[-1])
