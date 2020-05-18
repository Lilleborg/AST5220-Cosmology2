import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys

plt.style.use('bmh')

if "toy" in sys.argv[1:]:
    data_path = "../../data_toy/"
    figname = str(sys.argv[1]) + "_" + str(sys.argv[2])
else:
    data_path = "../../data/"
    figname = sys.argv[1]
filename = sys.argv[1] + ".txt"

quantities = np.loadtxt(data_path+filename,skiprows=1,max_rows=1,dtype=str)
data_toy = {}
data = {}
# loaded_data_toy = np.loadtxt(data_path_toy+filename,skiprows=1).T
loaded_data = np.loadtxt(data_path+filename,skiprows=2).T
for i,q in enumerate(quantities):
    # data_toy[q] = loaded_data_toy[i]
    data[q] = loaded_data[i]

# ell = data_toy["ell"]
# fig, ax = plt.subplots(figsize=(9,12))
# ax.plot(ell, data_toy["Cell"])
# ax.set_title("Toy cosmology")
# ax.set_yscale("log")
# ax.set_xscale("log")
# fig.savefig("../figs/CMB_powerspectrum_toy.png")
# fig.savefig("../figs/CMB_powerspectrum_toy.pdf")

ell = data["ell"]
fig, ax = plt.subplots(figsize=(9,12))
for i,q in enumerate(quantities[1:]):
    ax.plot(ell, data[q],label=q)
ax.set_title("LambdaCDM cosmology")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
fig.savefig("../figs/" +figname + ".png")

plt.show()