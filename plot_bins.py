import matplotlib.pyplot as plt
import numpy as np

# bin_vals = np.loadtxt("bins_output.txt")
#
# plt.bar(bin_vals[:,0], bin_vals[:,1])
# plt.show()
#
# gauss_vals = np.loadtxt("gauss_output.txt")
# plt.hist(gauss_vals)
# plt.show()

# pot_vals = np.loadtxt("pot_file.txt")
#
# for index in range(len(pot_vals[:, 1])):
#     if pot_vals[index, 1] > 1000000:
#         pot_vals[index, 1] = 0
#
#
# plt.plot(pot_vals[:, 0], pot_vals[:, 1])
# plt.show()

biased_vals = np.loadtxt("biased_pot_file.txt")
# plt.plot(biased_vals[:, 0], biased_vals[:, 1])
plt.plot(biased_vals[:, 0], biased_vals[:, 2], color="red")

plt.show()