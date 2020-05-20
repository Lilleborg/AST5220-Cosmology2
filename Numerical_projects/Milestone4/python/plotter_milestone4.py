import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys

plt.style.use('bmh')

if "toy" in sys.argv[1:]:
    data_path = "../../data_toy/"
    figname_C_ell = str(sys.argv[1]) + "_toy"
    figname_comp_PS = str(sys.argv[2]) + "_toy"
else:
    data_path = "../../data/"
    figname_C_ell = sys.argv[1]
    figname_comp_PS = sys.argv[2]
filename_C_ell = sys.argv[1] + ".txt"
filename_comp_PS = sys.argv[2] + ".txt"

# Cell stuff
quantities = np.loadtxt(data_path+filename_C_ell,skiprows=1,max_rows=1,dtype=str)
data = {}
loaded_data = np.loadtxt(data_path+filename_C_ell,skiprows=2).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i]

ell = data["ell"]
fig, ax = plt.subplots(figsize=(9,12))
for i,q in enumerate(quantities[1:]):
    ax.plot(ell, data[q],label=q)
ax.set_title("CMB power spectrum")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
# fig.savefig("../figs/" +figname_C_ell + ".png")

# Component power spectrum stuff
k_eq = np.loadtxt(data_path+filename_comp_PS,max_rows=1,usecols=1)
quantities = np.genfromtxt(data_path+filename_comp_PS,skip_header=1,max_rows=1,dtype=str, delimiter=",",autostrip=True)
data = {}
loaded_data = np.loadtxt(data_path+filename_comp_PS,skiprows=2).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i]

k = data[quantities[0]]
fig, ax = plt.subplots(2,1,figsize=(9,12),sharex=True)
for i,q in enumerate(quantities[1:]):
    if q == "radiation":
        ax[1].plot(k, data[q],label=q)
    else:
        ax[0].plot(k, data[q],label=q)
ax[0].axvline(x=k_eq,label=r"$k_{\rm{eq}} \approx$"+" {:.3f}".format(k_eq))
fig.suptitle("Component power spectrum")
ax[1].set_xlabel(r"$k [h/Mpc]$")
ax[0].set_ylabel(r"$P(k) [Mpc/h]^3$")
ax[1].set_ylabel(r"$P(k) [Mpc/h]^3$")
ax[0].set_yscale("log")
ax[1].set_yscale("log")
ax[0].set_xscale("log")
ax[1].set_xscale("log")
ax[0].legend()
ax[1].legend()
# fig.savefig("../figs/" +figname_comp_PS + ".png")

plt.show()