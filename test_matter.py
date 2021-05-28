import numpy as np
import matplotlib.pyplot as plt

matter = np.loadtxt('matter.txt')
k      = matter[:, 0]
P      = matter[:, 1]

nett   = np.loadtxt('matter_web.txt', skiprows=1)
k_nett = nett[:, 0]
P_nett = nett[:, 1]

data   = np.loadtxt('matter_data.txt', skiprows=1)
k_data = data[:, 0]
P_data = data[:, 1]
error  = data[:, 2]

# plt.loglog(k, P, label='Meg')
# plt.loglog(k_nett, P_nett, ls='--', label='Nett')
# plt.xlabel(r'$k$' + ' [h/Mpc]')
# plt.ylabel(r'$P(k)$' + ' [(Mpc/h)' + r'$^3$' + ']')
# plt.tight_layout()
# plt.legend()
# plt.show()

plt.plot(k, P, label='Meg')
plt.errorbar(k_data, P_data, error, fmt='x', label='data')
plt.xlabel(r'$k$' + ' [h/Mpc]')
plt.ylabel(r'$P(k)$' + ' [(Mpc/h)' + r'$^3$' + ']')
plt.tight_layout()
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()