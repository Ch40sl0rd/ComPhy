#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 11:01:18 2020

@author: lars
"""


import numpy as np
import matplotlib.pyplot as plt

file1 = 'Runge_Kutta4_diff.txt'
file2= 'ExterneKopplung.txt'
file3 = 'ExterneKraft'

t = np.linspace(0.0, 60.0, 10000)
diff = np.loadtxt(file1, usecols=0)
fig1 = plt.figure(1, figsize=(5,5))
ax1 = fig1.add_subplot(111)
ax1.plot(t, diff)
ax1.set_ylabel('Differenz zwischen RK4 und Exakter Lösung')
ax1.set_xlabel('Zeit t')
ax1.set_title('Runge-Kutta-Verfahren der Stufe 4')
ax1.grid()
fig1.savefig('Runge_Kutta_Proof.pdf')

diff = np.loadtxt(file2, usecols=0)
fig2 = plt.figure(2, figsize=(5,5))
ax2 = fig2.add_subplot(111)
ax2.plot(t, diff)
ax2.set_ylabel('Differenz der Lösungen')
ax2.set_xlabel('Zeit t')
ax2.grid()
ax2.set_title('Erzwungene Schwingung')
fig2.savefig('Erzwungene_Schwingung.pdf')

t = np.loadtxt(file3, usecols=0)
y = np.loadtxt(file3, usecols=1)
fig3 = plt.figure(3, figsize=(5,5))
ax3 = fig3.add_subplot(111)
ax3.plot(t, y)
ax3.set_ylabel('Position x(t)')
ax3.set_xlabel('Zeit t')
ax3.grid()
ax3.set_title('Externe Unbekannte Kraft')
fig3.savefig('unbekannteKraft.pdf')

n = np.arange(0,12)
avg_p = np.loadtxt('Power_Average_1.txt', usecols=0)
fig4 = plt.figure(4, figsize=(5,5))
ax4 = fig4.add_subplot(111)
ax4.plot(n, avg_p)
ax4.set_xlabel('Nummer der Periode')
ax4.set_ylabel('Durchschnittliche Leistung der Periode')
ax4.set_title('Mittlere Leistung im Einschwingvorgang')
ax4.grid()
fig4.savefig('Einschwing.pdf')

omega = np.loadtxt('Scan_Omega.txt', usecols=0)
P = np.loadtxt('Scan_Omega.txt', usecols=1)
fig5 = plt.figure(5, figsize=(5,5))
ax5 = fig5.add_subplot(111)
ax5.plot(omega, P)
ax5.set_xlabel('Eigenfrequenz $\omega_0$')
ax5.set_ylabel('Durchschnittliche Leistung')
ax5.set_title('Scan in $\omega_0$')
ax5.grid()
fig5.savefig('Scan_Omega_Power.pdf')

