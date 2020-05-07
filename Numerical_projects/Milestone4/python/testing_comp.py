import numpy as np
import matplotlib.pyplot as plt

data_path = "../../data_testing/"
k_values, _, horizon_entry = np.loadtxt(data_path+"perturbations_k_values.txt",skiprows=1,delimiter="|",unpack=True)
quantities = np.loadtxt(data_path+'component_test_k{:.5f}.txt'.format(k_values[0]),skiprows=1,max_rows=1,dtype=str)
data = {key: [] for key in quantities}
for k_value in k_values:
    loaded_data = np.loadtxt(data_path+'component_test_k{:.5f}.txt'.format(k_value),skiprows=2)
    for i,q in enumerate(quantities):
        data[q].append(loaded_data[:,i])

x = data["x"][0]

for ik in range(len(k_values)):
    fig,ax = plt.subplots()
    ax.set_title("k: "+str(k_values[ik]))
    ax.plot(x,data["dPsi"][ik],label="dPsi")
    ax.plot(x,data["dPhi"][ik],label="dPhi")
    ax.legend()
plt.show()