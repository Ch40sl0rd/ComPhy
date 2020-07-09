# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:49:17 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt

files = ['amp0_5.txt', 'amp1_0.txt', 'amp2_0.txt', 'amp4_0.txt']
amp = 0.25
for i in range(4):
    amp *=2
    t = np.loadtxt(files[i], usecols=0)
    e = np.loadtxt(files[i], usecols=1)
    plt.plot(t,e)
    plt.title('Verlauf des Echo-Signals mit $a=b={0}$'.format(amp))
    plt.xlabel('Zeit t')
    plt.ylabel('Signal $|e(t)|^2$')
    plt.grid()
    plt.savefig('aufgabe10_2_{0}.png'.format(i))
    plt.show()
    plt.close()
