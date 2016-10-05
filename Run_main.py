# Runs main.cpp and plots the generated .txt files

import os

os.system('c++ main.cpp jacobi.cpp non_interacting_case.cpp interacting_case.cpp unit_test.cpp -o main.o -O2 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')
os.system('./main.o')

from math import *
import numpy as np
from matplotlib import pyplot as plt

w_r = ['0.01', '0.5', '1', '5']
n = '350'

def read_rho_wavefunctions(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    rho = []; u1 = []; u2 = []; u3 = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        rho.append(float(words[0]))
        u1.append(float(words[1]))
        u2.append(float(words[2]))
        u3.append(float(words[3]))
    infile.close()
    return rho, u1, u2, u3

# Fetching data by a call on read_x_u_v for three different n:
rho, u1, u2, u3 = read_rho_wavefunctions('Results_non_interacting_case.txt')# + n[i] + '.txt')

# Plotting commands to look at the wave functions:
fig, ax = plt.subplots(1)
ax.set_title('Wave function with $n =$ ' + n)
ax.set_xlabel('$\\rho$')
ax.set_ylabel('$\mid u(\\rho) \mid^2$')
ax.plot(rho,u1,'r-',label='$\mid u_0(\\rho) \mid^2$')
ax.plot(rho,u2,'b-',label='$\mid u_1(\\rho) \mid^2$')
ax.plot(rho,u3,'g-',label='$\mid u_2(\\rho) \mid^2$')
ax.legend(loc='upper right',fancybox='True')
ax.grid()
plt.savefig('Non_interacting_case.pdf')

for i in range(0, 4):
    # Fetching data by a call on read_x_u_v for three different n:
    rho, u1, u2, u3 = read_rho_wavefunctions('Results_interacting_case_' + w_r[i] + '.txt')
    
    # Plotting commands to look at the wave functions:
    fig, ax = plt.subplots(1)
    ax.set_title('Wave function with $\omega_r =$ ' + w_r[i])
    ax.set_xlabel('$\\rho$')
    ax.set_ylabel('$\mid u(\\rho) \mid^2$')
    ax.plot(rho,u1,'r-',label='$\mid u_0(\\rho) \mid^2$')
    ax.plot(rho,u2,'b-',label='$\mid u_1(\\rho) \mid^2$')
    ax.plot(rho,u3,'g-',label='$\mid u_2(\\rho) \mid^2$')
    ax.legend(loc='upper right',fancybox='True')
    ax.grid()
    plt.savefig('Interacting_case_' + w_r[i] + '.pdf')
