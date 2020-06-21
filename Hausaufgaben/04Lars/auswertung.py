# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 14:32:12 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt

file1 = '9_3_variation_stützstellen.txt'
file2 = 'Lennard_jones_pot.txt'
file3 = 'variation_np.txt'
file4 = 'variation_ny.txt'
file5 = 'variation_pmax.txt'
file6 = 'lj_variation_np.txt'
file7 = 'lj_variation_ny.txt'
file8 = 'lj_variation_pmax.txt'
file9 = 'lj_variation_ymin.txt'

x = np.loadtxt(file1, usecols=0)
y = np.loadtxt(file1, usecols=1)
plt.plot(x,y, 'b.')
plt.title('Relativer Fehler der numerischen Lösung')
plt.xlabel('Anzahl der Stützstellen')
plt.ylabel('Relativer Fehler')
plt.grid()
plt.savefig('9_3_relative_abweichung.pdf')
plt.show()

x = np.loadtxt(file2, usecols=0)
y = np.loadtxt(file2, usecols=1)
plt.figure(figsize=(8,5))
plt.plot(x,y, 'b:')
plt.title('Lennard-Jones-Potential im Impulsraum')
plt.xlabel('Impuls p')
plt.ylabel(r'Potential $v^{l=0}(p,p)$')
plt.grid()
plt.savefig('lennard_jones_pot.pdf')
plt.show()

x1 = np.loadtxt(file3, usecols=0)
y1 = np.loadtxt(file3, usecols=1)
x2 = np.loadtxt(file4, usecols=0)
y2 = np.loadtxt(file4, usecols=1)
x3 = np.loadtxt(file5, usecols=0)
y3 = np.loadtxt(file5, usecols=1)
fig1 = plt.figure(1, figsize=(9,10))
ax1 = fig1.add_subplot(311)
ax1.grid()
ax1.plot(x1,y1)
ax1.set_xlabel('Anzahl der Impulsstützstellen $n_p$')
ax1.set_ylabel('relativer Fehler')
ax2 = fig1.add_subplot(312)
ax2.grid()
ax2.plot(x2,y2)
ax2.set_xlabel('Anzahl der Stützstellen $n_y$')
ax2.set_ylabel('relativer Fehler')
ax3 = fig1.add_subplot(313)
ax3.plot(x3, y3)
ax3.set_xlabel('Maximaler Impuls $p_{max}$')
ax3.grid()
ax3.set_ylabel('relativer Fehler')
fig1.tight_layout()
fig1.savefig('rel_error_params.pdf')
fig1.show()

x1 = np.loadtxt(file6, usecols=0)
y1 = np.loadtxt(file6, usecols=1)
x2 = np.loadtxt(file7, usecols=0)
y2 = np.loadtxt(file7, usecols=1)
x3 = np.loadtxt(file8, usecols=0)
y3 = np.loadtxt(file8, usecols=1)
x4 = np.loadtxt(file9, usecols=0)
y4 = np.loadtxt(file9, usecols=1)
fig2 = plt.figure(2, figsize=(9,12))
ax4 = fig2.add_subplot(411)
ax4.grid()
ax4.plot(x1,y1)
ax4.set_xlabel('Anzahl der Impulsstützstellen $n_p$')
ax4.set_ylabel('Energieeigenwert')
ax5 = fig2.add_subplot(412, sharey=ax4)
ax5.grid()
ax5.plot(x2,y2)
ax5.set_xlabel('Anzahl der Stützstellen $n_y$')
ax5.set_ylabel('Energieeigenwert')
ax6 = fig2.add_subplot(413, sharey=ax4)
ax6.plot(x3, y3)
ax6.set_xlabel('Maximaler Impuls $p_{max}$')
ax6.grid()
ax6.set_ylabel('Energieeignwert')

ax7 = fig2.add_subplot(414, sharey=ax4)
ax7.plot(x4,y4)
ax7.set_xlabel('Minimales $y_{min}$')
ax7.set_ylabel('Energieeigenwert')
ax7.grid()
fig2.tight_layout()
fig2.savefig('lj_variation_params.pdf')
fig2.show()