# -*- coding: utf-8 -*-
"""
Created on Fri May  1 14:34:59 2020
This script creates all tables and figures for the first homework
of Computerphysik SS2020.
@author: Lars Doepper
"""

import numpy as np
import matplotlib.pyplot as plt
import math

files = ['a_4.txt', 'a_2.txt', 'a_1.txt', 'a_0_1.txt', 'a_0_01.txt', 'a_0_001.txt']
fmts = ['r:', 'b:', 'g:', 'y:', 'k:', 'm:']
label = [r'a = 4', r'a = 2', r'a = 1', r'a = 0.1', r'a = 0.01', r'a = 0.001']

#files = ['a_4.txt', 'a_2.txt', 'a_1.txt',  'a_0_01.txt', 'a_0_001.txt']
#fmts = ['r:', 'b:', 'g:',  'k:', 'm:']
#label = [r'a = 4', r'a = 2', r'a = 1', r'a = 0.01', r'a = 0.001']

x = np.loadtxt(files[0], usecols=0)
y = np.empty((6,len(x)))

for i in range(6):
    y[i] = np.loadtxt(files[i], usecols=1)
    
fig = plt.figure(1, figsize=(5,4))
ax = fig.add_subplot()
for i in range(6):
    ax.plot(x, y[i], fmts[i], label=label[i])

ax.plot(x, 1/(math.sqrt(math.pi)*x), label='Referenzfunktion 1/z', color='burlywood')
ax.set_title('Elektrostatisches Potential')
ax.set_xlabel('Entfernung z')
ax.grid()
ax.set_ylabel('Potential V(z)')
ax.set_ylim(-0.5 , 7.5)
ax.legend()
fig.savefig('aufgabe1_referenzpot.pdf', dpi=400)


#create figure for electric field
e_field = np.loadtxt('E-Feld.txt', usecols=1)
e_field2 = np.loadtxt('E_Feld_calc.txt', usecols=1)
fig2 = plt.figure(2, figsize=(5,4))
ax2 = fig2.add_subplot()
ax2.plot(x, e_field,'b:', label='Elektrisches Feld aus numerischer Ableitung')
ax2.set_title('Elektrisches Feld')
ax2.set_xlabel('Entfernung z')
ax2.grid()
ax2.set_ylabel(r'Elektriches Feld $\vec{E}(z)$')
ax2.legend()
fig2.savefig('aufgabe2_1.pdf', dpi=400)
ax2.plot(x, e_field2,'g:', label='Elektrisches Feld aus numerischer Integration')
ax2.legend()
fig2.savefig('aufgabe2_2.pdf', dpi=400)

fig3 = plt.figure(3, figsize=(5,4))
z = np.loadtxt('e_field2_calc.txt', usecols=0)
e_z = np.loadtxt('e_field2_calc.txt', usecols=1)
ax3 = fig3.add_subplot()
ax3.plot(z, e_z, 'g:', label='Elektrisches Feld')
ax3.set_title('Elektrisches Feld von zwei Ladungsverteilungen')
ax3.set_xlabel('Entfernung z')
ax3.grid()
ax3.set_ylabel(r'Elektriches Feld $\vec{E}(z)$')
ax3.legend()
fig3.savefig('aufgabe3_1.pdf', dpi=400)

fig4 = plt.figure(4, figsize=(6,4))
a_2 = np.loadtxt('Nullstellen.txt', usecols=0)
nullstellen = np.loadtxt('Nullstellen.txt', usecols=1)
ax4 = fig4.add_subplot()
ax4.plot(a_2, nullstellen, 'b:', label="Gleichgewichtspunkte")
ax4.set_title('Gleichgewichtspunkte')
ax4.set_xlabel(r'Breite der Ladungsverteilung $a_2$')
ax4.set_ylabel('Nullstelle des Elektrischen Feldes')
ax4.grid()
ax4.legend()
fig4.savefig('Nullstellen.pdf', dpi=400)
