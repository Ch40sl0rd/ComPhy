# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:39:26 2020

@author: Lars
"""

import numpy as np
import matplotlib.pyplot as plt

file1 = "test.txt"

x = np.loadtxt(file1, usecols=0)
y = np.loadtxt(file1, usecols=1)

fig1 = plt.figure()
ax = fig1.add_subplot()
ax.plot(x,y)
