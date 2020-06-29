# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 10:35:27 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

def add_squares(x,y):
    return x*x+y*y

def solution(x,t):
    zaehler = 3+4*np.cosh(2*x-8*t)+np.cosh(4*x-64*t)
    nenner = (3*np.cosh(x-28*t)+np.cosh(3*x-36*t))
    return -12*zaehler/(nenner*nenner)

x = np.linspace(-30, 30, 10000)
t = np.linspace(-1, 1, 60)
X,T = np.meshgrid(x,t)
z = solution(X,T)

fig1 = plt.figure(num=1, figsize=(8,5))
axes1 = fig1.add_subplot(111)
cs = axes1.contourf(X,T,z, cmap='Blues')
fig1.colorbar(cs)
axes1.set_xlabel('Position x')
axes1.set_ylabel('Zeit t')
axes1.set_title('Contour-Plot der Solitonen')
fig1.savefig('aufgabe10_1_1.png')

fig2 = plt.figure(2, figsize=(7,5))
axes2 = fig2.add_subplot(111, projection='3d')
axes2.plot_surface(X,T,z)
axes2.set_xlabel('Position x')
axes2.set_ylabel('Zeit t')
#axes.set_zlabel('Abweichung von Normal Null')
axes2.view_init(70, 90)
axes2.set_title('Dreidimensionaler Verlauf der Solitonen')
fig2.savefig('aufgabe10_1_2.png')

fig3, axes3 = plt.subplots(figsize=(6,6))
for t in range(-1, 2):
    axes3.plot(x, solution(x,t), label='Zeit t={}'.format(t))
axes3.set_ylabel('Abweichung vom Normalpegel')
axes3.set_xlabel('Position x')
axes3.grid()
fig3.legend()
axes3.set_title('Zeitlicher Verlauf der Solitonen')
axes3.minorticks_on()
axes3.grid(which='major', linestyle='-', color='black')
axes3.grid(which='minor', linestyle=':', color='grey')
fig3.savefig('aufgabe10_1_3.png')
