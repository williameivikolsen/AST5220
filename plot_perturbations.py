import numpy as np
import matplotlib.pyplot as plt
import os
# import seaborn as sns

plt.rcParams.update({'font.size': 14})
# plt.style.use('seaborn')
# sns.set(font_scale=1.3)

data    = np.loadtxt('perturbations_k0.01.txt')

x       = data[:, 0]
theta0  = data[:, 1]
theta1  = data[:, 2]
theta2  = data[:, 3]
phi     = data[:, 4]
psi     = data[:, 5]

a = np.exp(x)

plt.plot(a, theta0, label=r'$\Theta_0$')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.show()

plt.plot(a, theta1, label=r'$\Theta_1$')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.show()

plt.plot(a, theta2, label=r'$\Theta_2$')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.show()

plt.plot(a, phi, label=r'$\Phi$')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.show()

plt.plot(a, psi, label=r'$\Psi$')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$x$')
plt.show()