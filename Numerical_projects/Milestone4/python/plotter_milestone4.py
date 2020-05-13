import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

plt.style.use('bmh')

data_path_toy = "../../data_toy/"
data_path = "../../data/"

quantities = np.loadtxt(data_path_toy+"Cells.txt",max_rows=1,dtype=str)
data_toy = {}
data = {}
loaded_data_toy = np.loadtxt(data_path_toy+"Cells.txt",skiprows=1).T
# loaded_data = np.loadtxt(data_path+"Cells.txt",skiprows=1).T
for i,q in enumerate(quantities):
    data_toy[q] = loaded_data_toy[i]
    # data[q] = loaded_data[i]

ell = data_toy["ell"]
fig, ax = plt.subplots(figsize=(9,12))
ax.plot(ell, data_toy["Cell"])
ax.set_title("Toy cosmology")
ax.set_yscale("log")
ax.set_xscale("log")
fig.savefig("../figs/CMB_powerspectrum_toy.png")
fig.savefig("../figs/CMB_powerspectrum_toy.pdf")

# fig, ax = plt.subplots(figsize=(9,12))
# ax.plot(ell, data["Cell"])
# ax.set_title("LambdaCDM cosmology")
# ax.set_yscale("log")
# ax.set_xscale("log")
# fig.savefig("../figs/CMB_powerspectrum.png")
# fig.savefig("../figs/CMB_powerspectrum.pdf")

plt.show()