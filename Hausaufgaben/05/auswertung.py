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