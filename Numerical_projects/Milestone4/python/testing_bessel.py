import numpy as np
import matplotlib.pyplot as plt

data_path = "../../data_testing/"
# Simple test of bessel function alone
j_ell, x = np.loadtxt(data_path+"testing_bessel.txt",unpack=True)

plt.plot(x,j_ell)

# Output from Milestone3 with source function
k_values, _, horizon_entry = np.loadtxt(data_path+"perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = ["x_array", "delta_cdm", "delta_b", "v_cdm", "v_b", "Phi", "Psi",\
    "Pi","Theta0", "Theta1", "Theta2","Source_T","arg","Source_T5","Source_T50","Source_T500","bessel50"]
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt(data_path+'perturbations_k{:.5f}.txt'.format(k_value),skiprows=1)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

x_array = data["x_array"][0]

id_ks = [0,1,2]
for ik in (id_ks):
    fig, ax = plt.subplots(3,1,figsize=(12,9),sharex=True)
    k = k_values[ik]
    ax[0].set_title("Integrand, k = {:.5f}".format(k))
    source_T5 = data["Source_T5"][ik]
    ax[0].plot(x_array,source_T5,label=r'$\ell = 5$')
    source_T50 = data["Source_T50"][ik]
    ax[0].plot(x_array,source_T50,label=r'$\ell = 50$')
    source_T500 = data["Source_T500"][ik]
    ax[0].plot(x_array,source_T500,label=r'$\ell = 500$')
    ax[0].legend()

    ax[1].set_title("Bessel alone, ell 50")
    ax[1].plot(x_array,data["bessel50"][ik])

    ax[2].set_title("Arg to bessel")
    ax[2].plot(x_array,data["arg"][ik])

plt.show()