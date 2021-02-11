import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn')

data = np.loadtxt('cosmology.txt')
x = data[:, 0]
eta = data[:, 1]
Hp = data[:, 2]
dHpdx = data[:, 3]
ddHpddx = data[:, 4]
OmegaB = data[:, 5]
OmegaCDM = data[:, 6]
OmegaLambda = data[:, 7]
OmegaR = data[:, 8]
OmegaNu = data[:, 9]
OmegaK = data[:, 10]
H = Hp/np.exp(x)

plt.plot(x, OmegaB + OmegaCDM, label=r'$\Omega_B + \Omega_{CDM}$')
plt.plot(x, OmegaLambda, label=r'$\Omega_{\Lambda}$')
plt.plot(x, OmegaR + OmegaNu, label=r'$\Omega_\gamma + \Omega_\nu$')
plt.plot(x, OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaNu + OmegaK, label=r'$\sum_i \Omega_i$')
plt.xlabel(r'$x$')
plt.legend()
plt.show()

plt.plot(x, H)
plt.xlabel(r'$x$')
plt.ylabel(r'$H(x)$')
plt.show()

plt.plot(x, Hp/H[0])
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathcal{H}(x)/H_0$')
plt.show()

plt.plot(x, dHpdx/Hp)
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$')
plt.show()

plt.plot(x, ddHpddx/Hp)
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$')
plt.show()

plt.plot(x, eta*Hp)
plt.xlabel(r'$x$')
plt.ylabel(r'$\eta(x)\mathcal{H}(x)$')
plt.show()

