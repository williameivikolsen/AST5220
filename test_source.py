import numpy as np
import matplotlib.pyplot as plt

data  = np.loadtxt('perturbations_k0.01.txt')
x     = data[:, 0]
S     = data[:, 10]
Sj5   = data[:, 11]
Sj50  = data[:, 12]
Sj500 = data[:, 13]
a     = np.exp(x)

plt.title(r'$S_k$' + ' for ' r'$k = 0.01$')
plt.plot(a, S)
plt.xscale('log')
plt.show()

plt.title(r'$S_k$' + ' for ' r'$k = 0.01$')
plt.plot(a, Sj5)
plt.xscale('log')
plt.show()

plt.title(r'$S_k$' + ' for ' r'$k = 0.01$')
plt.plot(a, Sj50)
plt.xscale('log')
plt.show()

plt.title(r'$S_k$' + ' for ' r'$k = 0.01$')
plt.plot(a, Sj500)
plt.xscale('log')
plt.show()