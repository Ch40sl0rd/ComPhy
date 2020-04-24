# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

name = 'ausgabe2.txt'
x = np.loadtxt(name, usecols=0)
y = np.loadtxt(name, usecols=1)

fig = plt.figure(1, figsize=(4,4))
ax1 = fig.add_subplot(111, label='Sinuskurve')
ax1.plot(x, y, 'b:')
#ax1.plot(x, np.sin(x), 'g:')
ax1.set_xlabel('x')
ax1.set_ylabel('sin(x)')
fig.show()