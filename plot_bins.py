import matplotlib.pyplot as plt
import numpy as np

bin_vals = np.loadtxt("bins_output.txt")

plt.bar(bin_vals[:,0], bin_vals[:,1])
plt.show()

gauss_vals = np.loadtxt("gauss_output.txt")
plt.hist(gauss_vals)
plt.show()