import numpy as np
import matplotlib.pyplot as plt

data_path = "../../data_testing/"
# # test of bessel function alone
# data_bessel = np.loadtxt(data_path+"testing_bessel.txt",skiprows=2).T
# labels = np.loadtxt(data_path+"testing_bessel.txt",skiprows=1,max_rows=1,dtype=str)
# x = data_bessel[0]
# arg = data_bessel[1]
# j_ell = [data_bessel[i] for i in range(2,len(data_bessel))]

# fig, ax = plt.subplots(2,1,figsize=(12,9),sharex=True)
# for i in range(len(j_ell)):
#     if i == 0:
#         ax[0].plot(x,arg,label="argument to bessel")
#     if np.all(j_ell[i]==0):
#         pass
#     elif i != 2:
#         pass
#     else:
#         ax[1].plot(x,j_ell[i],label=labels[i+2])
# ax[0].legend(loc=2)
# ax[1].legend(loc=2)

# Output from Milestone3 with source function
k_values, _, horizon_entry = np.loadtxt(data_path+"perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = np.loadtxt(data_path+'testing_perturbations_k{:.5f}.txt'.format(k_values[0]),skiprows=1,max_rows=1,dtype=str)
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt(data_path+'testing_perturbations_k{:.5f}.txt'.format(k_value),skiprows=2)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

x = data["x"][0]

for ik in range(len(k_values)):
    fig, ax = plt.subplots(4,1,figsize=(12,9),sharex=True)
    k = k_values[ik]
    fig.suptitle("k = "+str(k))
    
    ax[0].plot(x,data["ST"][ik],label="ST")
    ax[0].set_title(r"$\tilde{S}(k,x)$")
    ax[0].legend()

    for q in range(2,5):
        if np.all(data[quantities[q]][ik]==0):
            pass
        else:
            ax[1].plot(x,data[quantities[q]][ik]/1e-3,label=quantities[q])
    ax[1].set_title(r"$\tilde{S}(k,x)j_{\ell}[k(\eta_0-\eta(x]/10^{-3}$")
    ax[1].legend(loc=2)
    for ind_ell in range(10,len(quantities)):
        if np.all(data[quantities[ind_ell]][ik]==0) or ind_ell>13:
            pass
        else:
            ax[2].plot(x,data[quantities[ind_ell]][ik],label=quantities[ind_ell])
    ax[2].legend(loc=2)

    for q in range(6,10):
        ax[3].plot(x,data[quantities[q]][ik],label=quantities[q])
    ax[3].legend(loc=2)
plt.show()