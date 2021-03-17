import numpy as np
import matplotlib.pyplot as plt
import os
# import seaborn as sns

plt.rcParams.update({'font.size': 14})
# plt.style.use('seaborn')
# sns.set(font_scale=1.3)

h            = 0.67
OmegaB0      = 0.05
OmegaCDM0    = 0.267
TCMB0        = 2.7255
OmegaR0      = 5.50896e-05
OmegaLambda0 = 0.682945

data_background = np.loadtxt('cosmology.txt')
x = data_background [:, 0]
eta = data_background [:, 1]
Hp = data_background [:, 2]
dHpdx = data_background [:, 3]
ddHpddx = data_background [:, 4]
OmegaB = data_background [:, 5]
OmegaCDM = data_background [:, 6]
OmegaLambda = data_background [:, 7]
OmegaR = data_background [:, 8]
OmegaNu = data_background [:, 9]
OmegaK = data_background [:, 10]
H = Hp/np.exp(x)
a = np.exp(x)
c = 2.99792458e8
Mpc = 3.08567758e22

a_r = OmegaR0/(OmegaB0 + OmegaCDM0)                # Matter-radiation equality
a_d = ((OmegaB0 + OmegaCDM0)/OmegaLambda0)**(1/3)  # Matter-dark energy equality

supernova = np.loadtxt('supernova.txt', skiprows=1)
z_data = supernova[:, 0]
dL_data = supernova[:, 1]
error = supernova[:, 2]

data_recombination = np.loadtxt('recombination.txt')
x_recomb    = data_recombination[:, 0]
Xe          = data_recombination[:, 1]
log_ne      = data_recombination[:, 2]
tau         = data_recombination[:, 3]
dtaudx      = data_recombination[:, 4]
ddtauddx    = data_recombination[:, 5]
g           = data_recombination[:, 6]
dgdx        = data_recombination[:, 7]
ddgddx      = data_recombination[:, 8]

a_recomb = np.exp(x_recomb)

os.chdir('./Plots')

plt.figure()
plt.plot(a, OmegaB + OmegaCDM, label=r'$\Omega_B + \Omega_{CDM}$')
plt.plot(a, OmegaLambda, label=r'$\Omega_{\Lambda}$')
plt.plot(a, OmegaR + OmegaNu, label=r'$\Omega_\gamma + \Omega_\nu$')
plt.xlabel(r'$a$')
plt.ylabel(r'$\Omega$')
plt.xscale('log')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.legend()
plt.tight_layout()
# plt.savefig('omega.pdf')

plt.figure()
plt.loglog(a, H*Mpc/(100*1000))
plt.xlabel(r'$a$')
plt.ylabel(r'$H$' + ' ' + r'$\left[\frac{100 km/s}{Mpc}\right]$')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.legend()
plt.tight_layout()
# plt.savefig('H.pdf')

plt.figure()
plt.loglog(a, Hp*Mpc/(100*1000))
plt.xlabel(r'$a$')
plt.ylabel(r'$\mathcal{H}$' + ' ' + r'$\left[\frac{100 km/s}{Mpc}\right]$')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.legend()
plt.tight_layout()
# plt.savefig('Hp.pdf')

plt.figure()
plt.plot(a, dHpdx/Hp, label=r'$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$')
plt.plot(a, ddHpddx/Hp, label=r'$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$')
plt.xlabel(r'$a$')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.xscale('log')
plt.tight_layout()
plt.legend()
# plt.savefig('derivatives.pdf')

plt.figure()
plt.loglog(a, eta/Mpc)
plt.xlabel(r'$a$')
plt.ylabel(r'$\eta$' + ' [Mpc]')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.tight_layout()
plt.legend()
# plt.savefig('eta.pdf')

plt.figure()
plt.loglog(a, eta*Hp/c)
plt.xlabel(r'$a$')
plt.ylabel(r'$\eta\mathcal{H}/c$')
plt.axvline(x=1, ls='--', color='black', label=r'$a=1$')
plt.axvline(x=a_r, ls='--', color='red', label=r'$a_\gamma$')
plt.axvline(x=a_d, ls='--', color='magenta', label=r'$a_\Lambda$')
plt.legend()
plt.tight_layout()
# plt.savefig('etaHp_c.pdf')

eta0 = 4.39812e26
chi = (eta0 - eta)/(Mpc*1e3)
d_A = a*chi
d_L = d_A/a**2
z = 1/a - 1

plt.figure()
plt.plot(z[z <= np.max(z_data)], d_L[z <= np.max(z_data)], label=r'$d_L(z)$')
plt.errorbar(z_data, dL_data, error, fmt='.', markersize=9, label='Supernova data')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$z$')
plt.ylabel(r'$d_L$' + ' [Gpc]')
plt.legend()
plt.tight_layout()
# plt.savefig('luminosity_distance.pdf')

plt.figure()
plt.loglog(a_recomb, Xe)
plt.xlabel(r'$a$')
plt.ylabel(r'$X_e$')
plt.tight_layout()
plt.savefig('Xe.pdf')

plt.figure()
plt.loglog(a_recomb, np.exp(log_ne))
plt.xlabel(r'$a$')
plt.ylabel(r'$n_e$')
plt.tight_layout()
plt.savefig('n_e.pdf')

plt.figure()
plt.loglog(a_recomb, tau, label=r'$\tau$')
plt.loglog(a_recomb, -dtaudx, label=r'$-d\tau/dx$')
plt.loglog(a_recomb, ddtauddx, label=r'$d^2\tau/dx^2$')
plt.xlabel(r'$a$')
plt.tight_layout()
plt.legend()
plt.savefig('tau.pdf')

plt.figure()
plt.plot(a_recomb, g)
plt.xscale('log')
plt.xlabel(r'$a$')
plt.ylabel(r'$\tilde{g}$')
plt.tight_layout()
plt.savefig('g.pdf')

plt.figure()
plt.plot(a_recomb, dgdx)
plt.xscale('log')
plt.xlabel(r'$a$')
plt.ylabel(r'$d\tilde{g}/dx$')
plt.tight_layout()
plt.savefig('dgdx.pdf')

plt.figure()
plt.plot(a_recomb, ddgddx)
plt.xscale('log')
plt.xlabel(r'$a$')
plt.ylabel(r'$d^2\tilde{g}/dx^2$')
plt.tight_layout()
plt.savefig('ddgddx.pdf')

plt.show()