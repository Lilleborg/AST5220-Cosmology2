import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('./../../Milestone1/python')
from plotter_Milestone1 import color_each_regime

plt.style.use('bmh')

data_path = "../../data_testing/"

quantities = np.loadtxt(data_path+"test_bessel_spline.txt",max_rows=1,dtype=str)
data = {}
loaded_data = np.loadtxt(data_path+"test_bessel_spline.txt",skiprows=1).T
for i,q in enumerate(quantities):
    data[q] = loaded_data[i,::10]

x = data["arg"]
fig, ax = plt.subplots(figsize=(9,12))

for i,q in enumerate(quantities[len(quantities)-10:]):
    if i > 100:
        break
    ax.plot(x,data[q],label=r"$\ell = "+q+"$")
ax.legend()
plt.show()