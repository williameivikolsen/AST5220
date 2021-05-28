import numpy as np
import matplotlib.pyplot as plt
import os
# import seaborn as sns

plt.rcParams.update({'font.size': 14})
# plt.style.use('seaborn')
# sns.set(font_scale=1.3)

data = np.loadtxt('cells.txt')
ell  = data[:, 0]
cell = data[:, 1]

sw         = np.loadtxt('SW.txt')[:, 1]
isw        = np.loadtxt('ISW.txt')[:, 1]
doppler    = np.loadtxt('Doppler.txt')[:, 1]
quadrupole = np.loadtxt('Quadrupole.txt')[:, 1]

matter = np.loadtxt('matter.txt')
k      = matter[:, 0]
P      = matter[:, 1]

matter_data = np.loadtxt('matter_data.txt', skiprows=1)
k_data      = matter_data[:, 0]
P_data      = matter_data[:, 1]
P_error     = matter_data[:, 2]
k_eq        = 0.0201474 # h/Mpc

theta_data = np.loadtxt('theta.txt', skiprows=1)
k_theta    = theta_data[:, 0]
theta_100  = theta_data[:, 1]
theta_200  = theta_data[:, 2]
theta_500  = theta_data[:, 3]
theta_1000 = theta_data[:, 4]
eta0       = 14253.3 # Mpc

planck_low    = np.loadtxt('planck_cell_low.txt', skiprows=1)
ell_low       = planck_low[:, 0]
cell_low      = planck_low[:, 1]
err_up_low    = planck_low[:, 2]
err_down_low  = planck_low[:, 3]
planck_high   = np.loadtxt('planck_cell_high.txt', skiprows=1)
ell_high      = planck_high[:, 0]
cell_high     = planck_high[:, 1]
err_up_high   = planck_high[:, 2]
err_down_high = planck_high[:, 3]

os.chdir('./Plots')

plt.plot(k_theta, theta_100, label=r'$\ell = 100$')
plt.plot(k_theta, theta_200, label=r'$\ell = 200$')
plt.plot(k_theta, theta_500, label=r'$\ell = 500$')
plt.plot(k_theta, theta_1000, label=r'$\ell = 1000$')
plt.axvline(x=100/eta0, ls='--', color='black', label=r'$k = \ell/\eta_0$')
plt.axvline(x=200/eta0, ls='--', color='black')
plt.axvline(x=500/eta0, ls='--', color='black')
plt.axvline(x=1000/eta0, ls='--', color='black')
plt.xlabel(r'$k$' + ' [Mpc' + r'$^{-1}$' + ']')
plt.ylabel(r'$\Theta_\ell$')
plt.tight_layout()
plt.legend()
# plt.savefig('transferfunction.pdf')
plt.show()

plt.plot(k_theta, 100*(100+1)*theta_100**2/k_theta, label=r'$\ell = 100$')
plt.plot(k_theta, 200*(200+1)*theta_200**2/k_theta, label=r'$\ell = 200$')
plt.plot(k_theta, 500*(500+1)*theta_500**2/k_theta, label=r'$\ell = 500$')
plt.plot(k_theta, 1000*(1000+1)*theta_1000**2/k_theta, label=r'$\ell = 1000$')
plt.xlabel(r'$k$' + ' [Mpc' + r'$^{-1}$' + ']')
plt.ylabel(r'$\ell(\ell+1)\Theta_\ell^2/k$')
plt.tight_layout()
plt.legend()
# plt.savefig('integrand.pdf')
plt.show()

plt.loglog(ell, cell, label='Theory prediction')
plt.errorbar(ell_low, cell_low, yerr=[err_down_low, err_up_low], fmt='o', label='Planck 2018', color='r')
plt.errorbar(ell_high, cell_high, yerr=[err_down_high, err_up_high], fmt='o', color='r')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$' + ' [' + r'$\mu$' + 'K' + r'$^2$' + ']')
plt.tight_layout()
plt.legend()
# plt.savefig('powerspectrum.pdf')
plt.show()

plt.loglog(k, P, label='Theory prediction')
plt.axvline(x=k_eq, ls='--', label=r'$k_{eq}$', color='black')
plt.errorbar(k_data, P_data, P_error, fmt='o', label='SDSS Galaxies (DR7 LRG)', color='r')
plt.xlabel(r'$k$' + ' [h/Mpc]')
plt.ylabel(r'$P(k)$' + ' [(Mpc/h)' + r'$^3$' + ']')
plt.tight_layout()
plt.legend(frameon=True)
# plt.savefig('matterspectrum.pdf')
plt.show()

plt.loglog(ell, cell, label='Full S')
plt.loglog(ell, sw, label='SW')
plt.loglog(ell, isw, label='ISW')
plt.loglog(ell, doppler, label='Doppler')
plt.loglog(ell, quadrupole, label='Quadrupole')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$' + ' [' + r'$\mu$' + 'K' + r'$^2$' + ']')
plt.tight_layout()
plt.legend(frameon=True)
# plt.savefig('Source_func.pdf')
plt.show()