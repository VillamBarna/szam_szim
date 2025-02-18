import matplotlib.pyplot as plt
import numpy as np
import sho

sho.run_simulation(1, 0, 1, 1, 100, "asd", False)
data = np.loadtxt("energia_e.txt").T
print(data)
# plt.plot(data[0], 0.5*(np.power(data[1], 2)+np.power(data[2], 2)))
# plt.plot(data[0], data[1])
# plt.plot(data[1], data[2])
# plt.show()
