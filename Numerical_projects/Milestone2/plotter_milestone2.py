import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
sys.path.append('./../Milestone1/')
from plotter_Milestone1 import color_each_regime

# Plot styling
plt.style.use('bmh')
mpl.rc('font', family='serif', size=15)
mpl.rc('text', usetex=False)
figsize_single = (6,4.5)

# Load data
x_array, Xe, ne, tau, tau_deriv, tau_2deriv, g, g_deriv, g_2deriv =\
     np.loadtxt('../data/recombination.txt', unpack=True)

# Data handling and some numerical testing of g_tilde
x_array_ticks = np.append(np.linspace(x_array.min(), x_array.max(), 5),0)
integrated_g_tilde = np.trapz(g,x_array)
print('\nDifference between integrated g_tilde and 1:\nlog10(abs(integrated_g - 1)) =',\
     np.log10(np.abs(integrated_g_tilde-1)))
# g /= g_2deriv.max()
# g_deriv /=  g_2deriv.max()
# g_2deriv /= g_2deriv.max()    
g_deriv  /= (g_deriv.max()/g.max())/2    # Factor used to scale first derivative of g_tilde
g_2deriv /= (g_2deriv.max()/g.max())/4   # Factor used to scale second derivative of g_tilde
id_g_max = np.argmax(g)

all_axes = []
# Xe(x)
fig, Xeax = plt.subplots(figsize=figsize_single)
Xeax.plot(x_array,Xe)
Xeax.set_ylabel(r'$Xe = \frac{n_e}{n_b}$')
Xeax.set_title('Fraction of free electrons')
Xeax.axvline(x=-7.1649,linestyle='--',color='C4',label=r'$x_{rec}$',alpha=0.5)
Xeax.axhline(y=0.5,linestyle='-.',color='C5',label=r'$Xe = 0.5$',alpha=0.5)
Xeax.set_yscale('log')
all_axes.append(Xeax)

# tau(x)
fig, tauax = plt.subplots(figsize=figsize_single)
tauax.set_title('Optical depth and its derivatives')
tauax.plot(x_array,tau,label=r'$\tau$')
tauax.plot(x_array,-tau_deriv,label=r'$-\tau^{\prime}$')
tauax.plot(x_array,tau_2deriv,label=r'$\tau^{\prime\prime}$')
tauax.axvline(x=-6.98608,linestyle='--',color='C4',label=r'$x_{*}$',alpha=0.5)
tauax.set_yscale('log')
all_axes.append(tauax)

# g_tilde(x)
fig, gax = plt.subplots(2,figsize=(6,9))
fig.suptitle('Visibility function and\n its scaled derivatives')
for i in range(2):
    gax[i].plot(x_array,g,linewidth=3,alpha=0.7,label=r'$\tilde{g}$')
    gax[i].plot(x_array,g_deriv,linewidth=3,alpha=0.7,label=r'$\bar{\tilde{g}}^{\prime}$')
    gax[i].plot(x_array,g_2deriv,linewidth=3,alpha=0.7,label=r'$\bar{\tilde{g}}^{\prime\prime}$')
    gax[i].axvline(x=-6.98608,linestyle='--',color='C4',label=r'$x_{*}$',alpha=0.7)
    # gax[i].set_yscale('symlog')
all_axes.extend(gax)

for i,ax in enumerate(all_axes):
    color_each_regime(ax,x_array)
    ax.margins(x=0)
    ax.legend()
    ax.set_xticks(x_array_ticks)

# # gax[0].set_xlim(x_array[id_g_max]+x_array[id_g_max]*0.85,x_array[id_g_max]-x_array[id_g_max]*0.85)
gax[1].set_xlim(x_array[id_g_max]+x_array[id_g_max]*0.1,x_array[id_g_max]-x_array[id_g_max]*0.1)

# fig,axes = plt.subplots(2)
# axes[0].plot(x,g_deriv)
# axes[0].set_ylabel(r'$\tilde{g}^{\prime}$')
# axes[1].plot(x,g_2deriv)
# axes[1].set_ylabel(r'$\tilde{g}^{\prime\prime}$')



plt.show()