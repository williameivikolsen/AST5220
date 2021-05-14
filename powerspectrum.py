import numpy as np
import matplotlib.pyplot as plt
import os

data = np.loadtxt('cells.txt')
ell  = data[:, 0]
cell = data[:, 1]

theta = np.loadtxt('theta.txt')
k = 3000/0.67*theta[:, 0]
t_100 = theta[:, 1]
t_200 = theta[:, 2]
t_500 = theta[:, 3]
t_1000 = theta[:, 4]

os.chdir('./Plots')

plt.loglog(ell, cell)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
plt.savefig('powerspectrum.png')
plt.show()

plt.plot(k, t_100, label='100')
plt.plot(k, t_200, label='200')
plt.plot(k, t_500, label='500')
plt.plot(k, t_1000, label='1000')
plt.xlabel(r'$k$')
plt.ylabel(r'$\theta(k)$')
plt.legend()
plt.show()