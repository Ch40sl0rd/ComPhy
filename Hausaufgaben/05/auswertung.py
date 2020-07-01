# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 14:10:30 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt

def solution1(x,t):
    zaehler= -2.0
    nenner = np.cosh(x-4*t)**2
    return zaehler/nenner

def solution2(x,t):
    zaehler = 3+4*np.cosh(2*x-8*t)+np.cosh(4*x-64*t)
    nenner = (3*np.cosh(x-28*t)+np.cosh(3*x-36*t))
    return -12*zaehler/(nenner*nenner)

file1='funktion1.txt'
file2='funktion2.txt'
file3='funktion3.txt'

files=['funktion1.txt', 'funktion2.txt', 'funktion3.txt']

x = np.loadtxt(file1, usecols=0)
u = np.loadtxt(file1, usecols=1)
fig1 = plt.figure(1, figsize=(8,6))
ax1 = fig1.add_subplot(111)
ax1.plot(x,u, label='numerische Lösung')
ax1.plot(x, solution1(x,1.0), label='analytische Lösung')
ax1.set_title('Verlauf der Solitonen für N=1')
ax1.set_xlabel('Ort x')
ax1.set_ylabel('Abweichung von Standartpegel')
fig1.legend()
ax1.grid()
fig1.savefig('aufgabe10_2_1.png')
plt.show(fig1)
plt.close(fig1)

x = np.loadtxt(file2, usecols=0)
u = np.loadtxt(file2, usecols=1)
fig2 = plt.figure(2, figsize=(8,6))
ax2 = fig2.add_subplot(111)
ax2.plot(x,u, label='numerische Lösung')
ax2.plot(x, solution2(x,1.0), label='analytische Lösung')
ax2.set_title('Verlauf der Solitonen für N=2')
ax2.set_xlabel('Ort x')
ax2.set_ylabel('Abweichung von Standartpegel')
ax2.grid()
fig2.legend()
fig2.savefig('aufgabe10_2_2.png')
plt.show(fig2)
plt.close(fig2)

for i in range(3):
    x = np.loadtxt(files[i], usecols=0)
    u = np.loadtxt(files[i], usecols=1)
    plt.plot(x,u, label='numerische Lösung')
    plt.xlabel('Ort x')
    plt.ylabel('Abweichung von Normalpegel')
    plt.title('Lösung für N={0} und t=1'.format(i+1))
    plt.grid()
    plt.savefig('aufgabe10_N{0}.png'.format(i+1))
    plt.show()
    plt.close()

for i in range(1,10):
    x = np.loadtxt('funktion1_{0}.txt'.format(i), usecols=0)
    t = np.loadtxt('funktion1_{0}.txt'.format(i), usecols=1)
    plt.plot(x,t, label='N=1,{0}'.format(i))
plt.xlabel('Ort x')
plt.ylabel('Abweichung von Normalpegel')
plt.title('Nicht-ganzzahlige N bei t=1')
plt.legend()
plt.grid()
plt.savefig('aufgabe10_4.png')
plt.show()
plt.close()

x = np.loadtxt('ort.txt')
t = np.loadtxt('zeit.txt')
u = np.loadtxt('verlauf_funktion3.txt')
#Wir müssen hier den Background filtern, da sonst die
#erste Solitone nicht sichtbar ist.
np.putmask(u, np.abs(u)<1, 0.0)
np.putmask(u, np.abs(u)<3, u*5)
#Plotte den Contour-Plot der N=3 Solitonen
X,T = np.meshgrid(x,t)
plt.contourf(X,T,u, cmap='Blues')
plt.colorbar()
plt.xlabel('Ort x')
plt.ylabel('Zeit t')
plt.title('Verlauf der N=3 Lösung')
plt.savefig('aufgabe10_3_2.png')
plt.show()
plt.close()