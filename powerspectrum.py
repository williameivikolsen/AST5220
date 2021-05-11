import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('cells.txt')
ell  = data[:, 0]
cell = data[:, 1]

plt.loglog(ell, ell*(ell+1)*np.abs(cell)/(2*np.pi))
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)|C_\ell|/2\pi$')
plt.show()