import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('./../../Milestone1/python')
from plotter_Milestone1 import color_each_regime

plt.style.use('bmh')

data_path = "../../data/"

# Output from Milestone3 with source function
k_values, _, horizon_entry = np.loadtxt(data_path+"perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = np.loadtxt(data_path+'perturbations_k{:.5f}.txt'.format(k_values[0]),skiprows=1,max_rows=1,dtype=str)
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt(data_path+'perturbations_k{:.5f}.txt'.format(k_value),skiprows=2)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

x = data["x"][0]
# for ik in range(len(k_values)):
#     fig, ax = plt.subplots(4,1,figsize=(12,9),sharex=True)
#     k = k_values[ik]
#     fig.suptitle("k = "+str(k))
    
#     ax[0].plot(x,data["ST"][ik],label="ST")
#     ax[0].set_title(r"$\tilde{S}(k,x)$")
#     ax[0].legend()
#     color_each_regime(ax[0],x)

#     for q in range(2,5):
#         if np.all(data[quantities[q]][ik]==0):
#             pass
#         else:
#             ax[1].plot(x,data[quantities[q]][ik]/1e-3,label=quantities[q])
#     ax[1].set_title(r"$\tilde{S}(k,x)j_{\ell}[k(\eta_0-\eta(x]/10^{-3}$")
#     ax[1].legend(loc=2)
#     color_each_regime(ax[1],x)
#     for ind_ell in range(10,len(quantities)):
#         if np.all(data[quantities[ind_ell]][ik]==0) or ind_ell>13:
#             pass
#         else:
#             ax[2].plot(x,data[quantities[ind_ell]][ik],label=quantities[ind_ell])
#     ax[2].set_title("Bessel alone")
#     ax[2].legend(loc=2)
#     color_each_regime(ax[2],x)

#     for q in range(6,10):
#         ax[3].plot(x,data[quantities[q]][ik],label=quantities[q])
#     ax[3].set_title("Terms in source func")
#     ax[3].legend(loc=2)
#     color_each_regime(ax[3],x)
for ik in range(len(k_values)):
    fig, ax = plt.subplots()
    ax.plot(x,data["ST*j_5"][ik]/1e-3,label="ST*j_5")
    ax.plot(x,data["ST"][ik]*10,label="ST")
    # ax.plot(x,data["ST*j_50"][ik]/1e-3,label="ST*j_50")
    # ax.plot(x,data["ST*j_500"][ik]/1e-3,label="ST*j_500")
    
    ax.set_title("k: {:g}".format(k_values[ik]))
    ax.set_ylabel(r"$\tilde{S}(k,x)j_{\ell}[k(\eta_0-\eta(x]/10^{-3}$")
    ax.legend(loc=2)
    color_each_regime(ax,x)
    # fig.savefig('testing_perturbations_k{:.5f}.pdf'.format(k_values[ik]))
plt.show()