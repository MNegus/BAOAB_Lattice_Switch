import matplotlib.pyplot as plt
import numpy as np

bin_vals = np.loadtxt("bins_output.txt")

plt.bar(bin_vals[:,0], bin_vals[:,1])
plt.show()

gauss_vals = np.loadtxt("gauss_output.txt")
plt.hist(gauss_vals, bins=1000)
plt.show()



left_biased_vals = np.loadtxt("left_pot_file.txt")
plt.plot(left_biased_vals[:, 0], left_biased_vals[:, 2], color="red")

# plt.show()

right_biased_vals = np.loadtxt("right_pot_file.txt")
plt.plot(right_biased_vals[:, 0], right_biased_vals[:, 2])
plt.ylim(0, 1000)
plt.xlim(-10, 132)
plt.show()