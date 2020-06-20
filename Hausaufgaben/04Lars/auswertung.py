# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 14:32:12 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt

file1 = '9_3_variation_stützstellen.txt'
file2 = 'Lennard_jones_pot.txt'

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