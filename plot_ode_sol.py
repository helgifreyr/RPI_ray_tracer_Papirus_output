# Plot the solution that was generated

import matplotlib
matplotlib.use('Agg')
from numpy import loadtxt, pi
from pylab import figure, polar, plot, xlabel, ylabel, grid, hold, legend, title, savefig, show, subplot
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

tau, t, r, phi = loadtxt('output.dat', unpack=True)

figure(1, figsize=(6, 2.7))

# xlabel(r'$\phi$')
# ylabel(r'$r$')

# lw = 1


ax = subplot(111, polar=True)
ax.grid(False)
ax.axis('off')
ax.plot(phi, r,linewidth=0.8, color='black')
ax.set_ylim(0, 6)
initial = Circle((r[0], phi[0]), 0.02, transform=ax.transData._b, color='red', fill=True)
circle = Circle((0, 0), 1, transform=ax.transData._b, color='black', fill=True)
ax.add_artist(initial)
ax.add_artist(circle)

savefig('output.png', dpi=800)
