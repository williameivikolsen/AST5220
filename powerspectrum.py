import numpy as np
import matplotlib.pyplot as plt
import os
# import seaborn as sns

# plt.rcParams.update({'font.size': 14})
# plt.style.use('seaborn')
# sns.set(font_scale=1.3)

data = np.loadtxt('cells.txt')
ell  = data[:, 0]
cell = data[:, 1]

os.chdir('./test files')
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

os.chdir('..')
os.chdir('./Plots')

plt.loglog(ell, cell, label='Theory prediction')
plt.errorbar(ell_low, cell_low, yerr=[err_down_low, err_up_low], fmt='x', elinewidth=1, label='Planck data')
plt.errorbar(ell_high, cell_high, yerr=[err_down_high, err_up_high], fmt='x', elinewidth=1, label='Planck data')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
plt.tight_layout()
plt.legend()
plt.show()