import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys

plt.style.use('bmh')

if "toy" in sys.argv[1:]:
    data_path = "../../data_toy/"
    figname_C_ell = str(sys.argv[1]) + "_toy"
    figname_comp_PS = str(sys.argv[2]) + "_toy"
    figname_integrand = str(sys.argv[3]) + "_toy"
else:
    data_path = "../../data/"
    figname_C_ell = str(sys.argv[1])
    figname_comp_PS = str(sys.argv[2])
    figname_integrand = str(sys.argv[3])
filename_C_ell = sys.argv[1] + ".txt"
filename_comp_PS = sys.argv[2] + ".txt"
filename_integrand = sys.argv[3] + ".txt"

# Cell stuff
quantities = np.loadtxt(data_path+filename_C_ell,skiprows=1,max_rows=1,dtype=str)
data = {}
loaded_data = np.loadtxt(data_path+filename_C_ell,skiprows=2).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i]

ell = data["ell"]
fig, ax = plt.subplots(figsize=(6,4.5))
for i,q in enumerate(quantities[1:]):
    ax.plot(ell, data[q],label=q)
ax.set_title(r"CMB power spectrum, $C_\ell$")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel(r"$\ell$")
ax.set_ylabel(r"$C_\ell \, \frac{\ell(\ell+1)}{2\pi} \left(10^{6}\cdot T_{\rm{CMB,0}}\right)^2$")
ax.legend()
fig.savefig("../figs/" +figname_C_ell + ".pdf")

# Component power spectrum stuff
k_eq = np.loadtxt(data_path+filename_comp_PS,max_rows=1,usecols=1)
quantities = np.genfromtxt(data_path+filename_comp_PS,skip_header=1,max_rows=1,dtype=str, delimiter=",",autostrip=True)
data = {}
loaded_data = np.loadtxt(data_path+filename_comp_PS,skiprows=2).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i]

k = data[quantities[0]]
fig, ax = plt.subplots(2,1,figsize=(6,9),sharex=True)
for i,q in enumerate(quantities[1:]):
    if q == "radiation":
        ax[1].plot(k, data[q],label=q)
    else:
        ax[0].plot(k, data[q],label=q)
ax[0].axvline(x=k_eq,label=r"$k_{\rm{eq}} \approx$"+" {:.3f}".format(k_eq),linestyle=":")

ax[1].set_title("Radiation power spectrum")
ax[1].set_xlabel(r"$k, \, [h/Mpc]$")
ax[1].set_ylabel(r"$P(k), \, [Mpc/h]^3$")
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].legend()
ax[0].set_title("Matter power spectrum")
ax[0].set_ylabel(r"$P(k), \, [Mpc/h]^3$")
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[0].legend()
fig.savefig("../figs/" +figname_comp_PS + ".pdf")

# Integrand and transfer
ell_values = np.genfromtxt(data_path+filename_integrand,skip_header=1,max_rows=1,dtype=int)[2:]
quantities = np.genfromtxt(data_path+filename_integrand,skip_header=2,max_rows=1,dtype=str)
data = {}
loaded_data = np.loadtxt(data_path+filename_integrand,skiprows=3).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i]

k = data[quantities[0]]
Mpc = 3.08567758e22
k_Mpc = k*Mpc
fig, ax = plt.subplots(2,1,figsize=(6,9),sharex=True)
for i,q in enumerate(quantities[1:]):
    if "/k" not in q:
        ell = ell_values[i]
        ax[0].plot(k_Mpc,data[q],label=r"$\ell = $"+str(ell),alpha=0.7,linewidth=1)
        ax[1].plot(k_Mpc,ell*(ell+1)*data[q],label=r"$\ell = $"+str(ell),alpha=0.7,linewidth=1)

ax[0].set_ylabel(r"$\Theta_\ell(k)$")
ax[1].set_ylabel(r"$\ell(\ell+1)\Theta_\ell(k)$")
ax[1].set_xlabel(r"$k/Mpc$") 
handles, labels = ax[0].get_legend_handles_labels()
legend = fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
title = fig.suptitle("Transfer Functions")
print("Saving " + "../figs/"+figname_integrand+"_transfer.pdf")
fig.savefig("../figs/"+figname_integrand+"_transfer.pdf",bbox_extra_artists=(legend,title),bbox_inches='tight')

fig, ax = plt.subplots(2,1,figsize=(6,9),sharex=True)
for i,q in enumerate(quantities[1:]):
    if "/k" in q:
        ell = ell_values[i-len(ell_values)]
        ax[0].plot(k_Mpc,data[q]/Mpc,label=r"$\ell = $"+str(ell),alpha=0.7,linewidth=1)
        ax[1].plot(k_Mpc,ell*(ell+1)*data[q]/Mpc,label=r"$\ell = $"+str(ell),alpha=0.7,linewidth=1)

ax[0].set_ylabel(r"$\Theta_\ell(k)^2/k$")
ax[1].set_ylabel(r"$\ell(\ell+1)\Theta_\ell(k)^2/k$")
ax[1].set_xlabel(r"$k/Mpc$")
handles, labels = ax[0].get_legend_handles_labels()
legend = fig.legend(handles,labels,bbox_to_anchor=(1.0, 0.5), loc='center left')
title = fig.suptitle("Integrand")
print("Saving " + "../figs/"+figname_integrand+"_integrand.pdf")
fig.savefig("../figs/"+figname_integrand+"_integrand.pdf",bbox_extra_artists=(legend,title),bbox_inches='tight')

plt.show()