import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
sys.path.append('./../../Milestone1/python')
from plotter_Milestone1 import color_each_regime

# Plot styling
plt.style.use('bmh')
mpl.rc('font', family='serif', size=15)
mpl.rc('text', usetex=False)

# Load data
k_values, _, horizon_entry = np.loadtxt("../../data/perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = ["x_array", "delta_cdm", "delta_b", "v_cdm", "v_b", "Phi", "Psi",\
    "Pi","Theta0", "Theta1", "Theta2","Source_T","Source_T5","Source_T50","Source_T500"]
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt('../../data/perturbations_k{:.5f}.txt'.format(k_value),skiprows=1)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

    
x_array = data["x_array"][0]
x_array_ticks = np.append(np.linspace(x_array.min(), x_array.max(), 6), 0)
x_array_ticks = np.array(x_array_ticks,dtype=int)
all_axes = []

# Matter perturbations
matter_fig, matter_axes = plt.subplots(2,1,sharex=True,figsize=(6,9))
delta_ax, v_ax = matter_axes
for ik,k in enumerate(k_values):
    delta_ax.plot(x_array,np.abs(data["delta_cdm"][ik]),color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$',linewidth=1.3,alpha=0.9)
    delta_ax.plot(x_array,np.abs(data["delta_b"][ik]),'--',color='C{:0d}'.format(ik),linewidth=1.3,alpha=0.9)

    v_ax.plot(x_array,np.abs(data["v_cdm"][ik]),color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$',linewidth=1.3,alpha=0.9)
    v_ax.plot(x_array,np.abs(data["v_b"][ik]),'--',color='C{:0d}'.format(ik),linewidth=1.3,alpha=0.9)

handles, labels = delta_ax.get_legend_handles_labels()
matter_legend = matter_fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
matter_title = matter_fig.suptitle("Matter perturbations")
delta_ax.set_ylabel(r'$|\delta_{\rm{cdm}}|$, $|\delta_{\rm{b}}|$ in dashed')
delta_ax.set_yscale('log')
v_ax.set_ylabel(r'$|v_{\rm{cdm}}|$, $|v_{\rm{b}}|$ in dashed')
v_ax.set_yscale('log')
v_ax.set_xlabel(r'$x=\ln(a)$')
all_axes.extend(matter_axes)

# Multipoles
theta_fig, theta_axes = plt.subplots(2,1,sharex=True,figsize=(6,9))
theta0_ax,theta1_ax = theta_axes
for ik,k in enumerate(k_values):
    theta0_ax.plot(x_array,4*data["Theta0"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$',linewidth=1.3,alpha=0.9)
    theta1_ax.plot(x_array,-3*data["Theta1"][ik],linewidth=1.3,alpha=0.9)#,label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')

handles, labels = theta0_ax.get_legend_handles_labels()
theta_legend = theta_fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
theta_title = theta_fig.suptitle("Multipole perturbations,\n "+r"$\ell = 0,\,\ell=1$")
theta0_ax.set_ylabel(r'$\delta_{\gamma} = 4\theta_0$')
theta1_ax.set_ylabel(r'$v_{\gamma} = -3\theta_1$')
all_axes.extend(theta_axes)

# Potentials
potential_fig,potential_ax = plt.subplots(figsize=(6,4.5))
for ik,k in enumerate(k_values):
    potential_ax.plot(x_array,data["Phi"][ik],color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$',linewidth=1.3,alpha=0.9)
    potential_ax.plot(x_array,data["Psi"][ik],'--',color='C{:0d}'.format(ik),linewidth=1.3,alpha=0.9)
handles, labels = potential_ax.get_legend_handles_labels()
potential_legend = potential_fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
potential_title = potential_fig.suptitle("Gravitational potentials")
potential_ax.set_ylabel(r"$\Phi,\, \Psi$ in dashed")
all_axes.append(potential_ax)

# # Source function, not used in this milestone
SourceT_fig, SourceT_ax = plt.subplots()
SourceT_fig.suptitle(r'Temperature Source function')
for ik,k in enumerate(k_values):
    SourceT_ax.plot(x_array,data["Source_T"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
all_axes.append(SourceT_ax)

# Line of sight integrand
los_fig, los_ax = plt.subplots()
for ik, k in enumerate(k_values):
    los_ax.plot(x_array[x_array<0],data["Source_T5"][ik][x_array<0])
all_axes.append(los_ax)

for ax in all_axes:
    color_each_regime(ax,x_array)
    ax.margins(x=0)
    ax.set_xticks(x_array_ticks)
    for ik in range(len(k_values)):
        ax.axvline(x=horizon_entry[ik],linestyle=':',color='C{:0d}'.format(ik))

# Saving
print("Saving ../figs/matter_pert.pdf")
matter_fig.savefig("../figs/matter_pert.pdf",bbox_extra_artists=(matter_legend,matter_title),bbox_inches='tight')
print("Saving ../figs/multipole_pert.pdf")
theta_fig.savefig("../figs/multipole_pert.pdf",bbox_extra_artists=(theta_legend,theta_title),bbox_inches='tight')
print("Saving ../figs/potential_pert.pdf")
potential_fig.savefig("../figs/potential_pert.pdf",bbox_extra_artists=(potential_legend,potential_title),bbox_inches='tight')

plt.show()
