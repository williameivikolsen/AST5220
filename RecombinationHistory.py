import numpy as np
import matplotlib.pyplot as plt
#Sets the fontsize of plots, also use linewidth 5
#plt.rcParams.update({'font.size': 48})

#[0]=x, [1]=Xe, [2]=ne, [3]=tau, [4]=dtau/dx, [5]=ddtau/dxdx, [6]=g_tilde, [7]=dg_tilde/dx, [8]=ddg_tilde/dxdx,
data = np.loadtxt('recombination.txt')

'''Plots of Xe'''
figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 1], label = r"$X_e$", linewidth = 1)
axes.set_title (r"Evolution of the electron fraction")
axes.set_ylim (bottom = 10**(-4), top = 10**1)
axes.set_xlim (left = -12, right = 0)
axes.set_yscale ("log")
axes.set_xlabel (r"x")
axes.legend()
plt.show ()


'''Plots of g and its derivatives'''
figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 6], label = r"$\tilde{g}$", linewidth = 1)
axes.set_title (r"Evolution of the visibility function")
axes.set_ylim (bottom = -0.05, top = 5.05)
axes.set_xlim (left = -12, right = 0)
axes.set_xlabel (r"x")
axes.legend()
plt.show ()

figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 7], label = r"$\tilde{g}$", linewidth = 1)
axes.set_title (r"Evolution of the visibility function")
axes.set_ylim (bottom = min(data[:, 7]*1.05), top = max(data[:, 7])*1.05)
axes.set_xlim (left = -12, right = 0)
axes.set_xlabel (r"x")
axes.legend()
plt.show ()

figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 8], label = r"$\tilde{g}$", linewidth = 1)
axes.set_title (r"Evolution of the visibility function")
axes.set_ylim (bottom = min(data[:, 8])*1.01, top = max(data[:, 8])*1.01)
axes.set_xlim (left = -12, right = 0)
axes.set_xlabel (r"x")
axes.legend()
plt.show ()


'''Plots of tau and its derivatives'''
figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 3], label = r"$\tau$", linewidth = 1)
axes.plot (data[:, 0], -1*data[:, 4], label = r"$-\frac{d\tau}{dx}$", linewidth = 1)
axes.plot (data[:, 0], data[:, 5], label = r"$\frac{d^2\tau}{dx^2}$", linewidth = 1)
axes.set_title (r"Evolution of the optical depth")
axes.set_ylim (bottom = 10**(-8), top = 10**8)
axes.set_xlim (left = -12, right = 0)
axes.set_yscale ("log")
axes.set_xlabel (r"x")
axes.legend()
plt.show ()

# Plot of ne
figure, axes = plt.subplots ()
figure.set_size_inches (8, 8)
axes.plot (data[:, 0], data[:, 2], label = r"$n_e$", linewidth = 1)
axes.set_title (r"Evolution of electron number")
axes.set_xlabel (r"x")
axes.legend()
plt.show()