import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('./../../Milestone1/python')
from plotter_Milestone1 import color_each_regime

plt.style.use('bmh')

data_path = "../../data_testing/"
k_values, _, horizon_entry = np.loadtxt(data_path+"perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = np.loadtxt(data_path+'component_test_k{:.5f}.txt'.format(k_values[0]),skiprows=1,max_rows=1,dtype=str)
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt(data_path+'component_test_k{:.5f}.txt'.format(k_value),skiprows=2)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

x = data["x"][0]

figs = []
# for ik in range(len(k_values)):
#     fig,ax = plt.subplots()
#     ax.set_title("k: "+str(k_values[ik]))
#     ax.plot(x,data["dPsi_dx"][ik],label="dPsi_dx")
#     ax.plot(x,data["dPhi_dx"][ik],label="dPhi_dx")
#     ax.legend()
#     color_each_regime(ax,x)
#     figs.append(fig)

for ik in range(len(k_values)):
    fig,ax = plt.subplots()
    ax.set_title("k: "+str(k_values[ik]))
    ax.plot(x,data["term3"][ik],label="term3")
    ax.plot(x,data["term3_1"][ik],label="term3_1")
    ax.plot(x,data["term3_2"][ik],label="term3_2")
    ax.plot(x,data["term3_3"][ik],label="term3_3")
    ax.legend()
    color_each_regime(ax,x)
    figs.append(fig)
plt.show()
for ik,fig in enumerate(figs):
    fig.savefig("comp_test_k_{:.5f}.pdf".format(k_values[ik]))