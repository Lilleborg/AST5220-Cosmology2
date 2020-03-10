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

# Load data
x_array, Xe, Xe_saha, ne, tau, tau_deriv, tau_2deriv, g, g_deriv, g_2deriv =\
     np.loadtxt('../data/recombination.txt', unpack=True)

# Data handling and some numerical testing of g_tilde
# xstar and xrec taken from print from ./cmb
xstar      = -6.98608
xrec       = -7.1649
xrec_saha  = -7.23022

zoomed_xlim = [xstar+xstar*0.1, xstar-xstar*0.1]
x_array_ticks = np.append(np.linspace(x_array.min(), x_array.max(), 5),0)
x_zoomed_ticks = zoomed_xlim
x_zoomed_ticks.extend([xrec,xstar])

# Cut of the Saha approximation at some small value
Xe_saha = Xe_saha[Xe_saha>1e-5]
x_array_saha = x_array[:len(Xe_saha)]

# Check that unscaled g_tilde is a PDF
integrated_g_tilde = np.trapz(g,x_array)
print('\nDifference between integrated g_tilde and 1:\nlog10(abs(integrated_g - 1)) =',\
     np.log10(np.abs(integrated_g_tilde-1)))

# Factor to scale the square of absolute value of the visibility functions to be PDFs, like the wavefunction in QM
scalefactor_g        = (np.trapz(np.abs(g)**2       ,x_array))**(1/2)
scalefactor_g_deriv  = (np.trapz(np.abs(g_deriv)**2 ,x_array))**(1/2)
scalefactor_g_2deriv = (np.trapz(np.abs(g_2deriv)**2,x_array))**(1/2)

print('\nCheck that the squared absolute value of the visibility functions are PDFs:')
print('integrated |g|^2:       ',np.trapz(np.abs(g/scalefactor_g)**2,x_array))
print('integrated |g_deriv|^2: ',np.trapz(np.abs(g_deriv/scalefactor_g_deriv)**2,x_array))
print('integrated |g_2deriv|^2:',np.trapz(np.abs(g_2deriv/scalefactor_g_2deriv)**2,x_array))

all_axes = []
# Xe(x)
Xefig, Xeax = plt.subplots(1,2,figsize=(12,4.5))
all_axes.extend(Xeax)
for i in range(2):
    Xeax[i].semilogy(x_array,Xe,label='$Xe$')
    Xeax[i].semilogy(x_array_saha,Xe_saha,ls='-.',label='$Xe_{Saha}$',color='C5')
    Xeax[i].tick_params(axis='y',labelcolor='C0')
    Xeax[i].axhline(y=0.5,label=r'$Xe=0.5$',linestyle='-.',color='C0',linewidth=1)
    Xeax[i].set_ylim(1e-4,1.4)
    
    # ne(x) overplotted
    neax = Xeax[i].twinx()
    neplot = neax.semilogy(x_array,ne,color='C1',linewidth=1)
    neax.tick_params(axis='y',labelcolor='C1')

    if i == 0:
        Xeax[i].set_ylabel(r'$Xe = \frac{n_e}{n_b}$',color='C0')
        Xeax[i].set_xticks(x_array_ticks)
        Xeax[i].set_xticklabels(x_array_ticks,rotation=10)
    else:
        neax.set_ylabel(r'$n_e \, [\rm{m^{-3}}] $',color='C1')
        Xeax[i].set_xticks(x_zoomed_ticks)
        Xeax[i].set_xticklabels(['{:.2f}'.format(i) for i in x_zoomed_ticks],rotation=10)

# tau(x)
taufig, tauaxes = plt.subplots(1,2,figsize=(12,4.5))
all_axes.extend(tauaxes)
for i in range(2):
    tauaxes[i].semilogy(x_array,tau,label=r'$\tau$')
    tauaxes[i].semilogy(x_array,-tau_deriv,label=r'$-\tau^{\prime}$')
    tauaxes[i].semilogy(x_array,tau_2deriv,label=r'$\tau^{\prime\prime}$')

    # unscaled visibility function overplotted
    g_tau_ax = tauaxes[i].twinx()
    g_tau_plot = g_tau_ax.plot(x_array,g,color='C3',alpha=0.7,linewidth=1)
    g_tau_ax.tick_params(axis='y',labelcolor='C3')

    if i == 0:
        tauaxes[i].set_xticks(x_array_ticks)
        tauaxes[i].set_xticklabels(x_array_ticks,rotation=10)
        tauaxes[i].set_ylabel(r'$\tau$')
    else:
        g_tau_ax.set_ylabel(r'$\tilde{g}$',color='C3')
        tauaxes[i].set_xticks(x_zoomed_ticks)
        tauaxes[i].set_xticklabels(['{:.2f}'.format(i) for i in x_zoomed_ticks],rotation=10)

# g_tilde(x)
gfig, gax = plt.subplots(1,2,figsize=(12,4.5))
all_axes.extend(gax)
for i in range(2):
    gax[i].plot(x_array,g/scalefactor_g,alpha=0.7,label=r'$\overline{\tilde{g}}$')
    gax[i].plot(x_array,g_deriv/scalefactor_g_deriv,alpha=0.7,label=r'$\overline{\tilde{g}^{\prime}}$')
    gax[i].plot(x_array,g_2deriv/scalefactor_g_2deriv,alpha=0.7,label=r'$\overline{\tilde{g}^{\prime\prime}}$')
    
    if i == 0:
        gax[i].set_xticks(x_array_ticks)
        gax[i].set_xticklabels(x_array_ticks,rotation=10)
    else:
        gax[i].set_xticks(x_zoomed_ticks)
        gax[i].set_xticklabels(['{:.2f}'.format(i) for i in x_zoomed_ticks],rotation=10)

# Some common handling of all axes
for i,ax in enumerate(all_axes):
    color_each_regime(ax,x_array)
    ax.axvline(x=xstar,linestyle='--',color='C3',label=r'$x_{*}$',linewidth=1)
    ax.axvline(x=xrec,linestyle=':',color='C4',label=r'$x_{rec}$',linewidth=1)
    ax.margins(x=0)

# Final tweaking and saving
# Xe and ne
Xeax[1].set_xlim(*zoomed_xlim)
handles, labels = Xeax[0].get_legend_handles_labels()
Xelegend = Xefig.legend(handles+neplot,labels+['\n'+'$n_e$'],bbox_to_anchor=(1.0,0.5),loc='center left')
Xetitle = Xefig.suptitle('Fraction of free electrons and the number density')
Xefig.savefig('./figs/free_electrons.pdf',bbox_extra_artists=(Xelegend,Xetitle),bbox_inches='tight')

# Tau
tauaxes[1].set_xlim(*zoomed_xlim)
handles, labels = tauaxes[0].get_legend_handles_labels()
taulegend = taufig.legend(handles+g_tau_plot,labels+['\n'+r'$\tilde{g}$'],bbox_to_anchor=(1.0,0.5),loc='center left')
tautitle = taufig.suptitle('Optical depth, its derivatives and the visibility function')
taufig.savefig('./figs/optical_depth.pdf',bbox_extra_artists=(taulegend,tautitle),bbox_inches='tight')

# g_tilde
gax[1].set_xlim(*zoomed_xlim)
handles, labels = gax[0].get_legend_handles_labels()
glegend = gfig.legend(handles,labels,bbox_to_anchor=(1.0,0.5),loc='center left')
gtitle = gfig.suptitle('Scaled visibility function and derivatives')
gfig.savefig('./figs/visibility_functions.pdf',bbox_extra_artists=(glegend,gtitle),bbox_inches='tight')
