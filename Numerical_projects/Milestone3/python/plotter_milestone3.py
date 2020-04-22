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
    "Pi","Theta0", "Theta1", "Theta2","Source_T"]
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
    delta_ax.plot(x_array,np.abs(data["delta_cdm"][ik]),color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
    delta_ax.plot(x_array,np.abs(data["delta_b"][ik]),'--',color='C{:0d}'.format(ik))

    v_ax.plot(x_array,np.abs(data["v_cdm"][ik]),color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
    v_ax.plot(x_array,np.abs(data["v_b"][ik]),'--',color='C{:0d}'.format(ik))

handles, labels = delta_ax.get_legend_handles_labels()
matter_legend = matter_fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
matter_title = matter_fig.suptitle("Matter perturbations")
delta_ax.set_title(r'$|\delta_{\rm{cdm}}|$, $|\delta_{\rm{b}}|$ in dashed')
delta_ax.set_yscale('log')
v_ax.set_title(r'$|v_{\rm{cdm}}|$, $|v_{\rm{b}}|$ in dashed')
v_ax.set_yscale('log')
v_ax.set_xlabel(r'$x=\ln(a)$')
all_axes.extend(matter_axes)

# Multipoles
theta_fig, theta_axes = plt.subplots(2,1,sharex=True,figsize=(6,9))
theta0_ax,theta1_ax = theta_axes
for ik,k in enumerate(k_values):
    theta0_ax.plot(x_array,4*data["Theta0"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
    theta1_ax.plot(x_array,-3*data["Theta1"][ik])#,label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')

handles, labels = theta0_ax.get_legend_handles_labels()
theta_legend = theta_fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
theta_title = theta_fig.suptitle("Multipole perturbations, "+r"$\ell = 0,\,\ell=1$")
theta0_ax.set_title(r'$\delta_{\gamma} = 4\theta_0$')
theta1_ax.set_title(r'$v_{\gamma} = -3\theta_1$')
all_axes.extend(theta_axes)
# potential_fig,potential_ax = plt.subplots()
# potential_fig.suptitle(r'$\Phi$, $\Psi$ in dashed')
# for ik,k in enumerate(k_values):
#     potential_ax.plot(x_array,data["Phi"][ik],color='C{:0d}'.format(ik),label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
#     potential_ax.plot(x_array,data["Psi"][ik],'--',color='C{:0d}'.format(ik))
# all_axes.append(potential_ax)

# theta0_fig, theta0_ax = plt.subplots()
# theta0_fig.suptitle(r'$\delta_{\gamma} = 4\theta_0$')
# for ik,k in enumerate(k_values):
#     theta0_ax.plot(x_array,4*data["Theta0"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
# all_axes.append(theta0_ax)

# theta1_fig, theta1_ax = plt.subplots()
# theta1_fig.suptitle(r'$v_{\gamma} = -3\theta_1$')
# for ik,k in enumerate(k_values):
#     theta1_ax.plot(x_array,-3*data["Theta1"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
# all_axes.append(theta1_ax)

# SourceT_fig, SourceT_ax = plt.subplots()
# SourceT_fig.suptitle(r'Temperature Source function')
# for ik,k in enumerate(k_values):
#     SourceT_ax.plot(x_array,data["Source_T"][ik],label=r'$k = {:g}$'.format(k)+r'$/\rm{Mpc}$')
# all_axes.append(SourceT_ax)

for ax in all_axes:
    color_each_regime(ax,x_array)
    ax.margins(x=0)
    # ax.legend()
    ax.set_xticks(x_array_ticks)
    for ik in range(len(k_values)):
        ax.axvline(x=horizon_entry[ik],linestyle=':',color='C{:0d}'.format(ik))
plt.show()