import matplotlib.pyplot as plt
import numpy as np

vx, vy, vz = np.loadtxt('Velocities.txt',delimiter=' ', unpack=True)

fig = plt.figure(1)
ax1 = fig.add_subplot(311)
plt.subplots_adjust(hspace=0.4)
plt.title("$v_x$")
ax1.hist(vx, bins=30)

ax2 = fig.add_subplot(312)
plt.title("$v_y$")
ax2.hist(vy, bins=30)

ax3 = fig.add_subplot(313)
plt.title("$v_z$")
ax3.hist(vz, bins=30)

plt.show()