#!/usr/bin/python
# -*- coding: utf-8 -*- 
import random
import math
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import random
from scipy.linalg import expm
from numpy import linalg as LA
from numpy.linalg import matrix_power
import scipy as sp
import scipy.linalg
import time 
import datetime 
import sys
startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S"))

if len(sys.argv) < 5:
  print("Usage:", str(sys.argv[0]), "READIN " " SAVE_or_NOT " "NCOL " "ITERS")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])
NCOL = int(sys.argv[3])
Niters_sim = int(sys.argv[4])
NMAT = 2
g=1.
h=1.
GAP=1
dt=1e-4
nsteps = int(1e-2/dt)
cut=int(0.25*Niters_sim) 
#•••••••••••••••••••••••••••••••••••••••••••••••••

X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM = []
expDS = [] 
EOS = [] 
trX2 = [] 
trX4 = []
MOM = []
monitor = []
evals = [] 

print ("Hoppe 2-matrix model")
print ("NCOL=" "%3.0f " ","  " and '(g = h = )' = " " (%4.2f)" % (NCOL, g))
print ("---------------------------------------------------------------------------------")

def dagger(a):
    return np.transpose(a).conj()

def box_muller():
    PI = 2.0*math.asin(1.0);
    r = random.uniform(0,1)
    s = random.uniform(0,1)
    p = np.sqrt(-2.0*np.log(r)) * math.sin(2.0*PI*s)
    q = np.sqrt(-2.0*np.log(r)) * math.cos(2.0*PI*s)
    return p,q

def comm(A,B):
    return np.dot(A,B) - np.dot(B,A)

def unit_matrix():
    matrix = np.zeros((NCOL, NCOL), dtype=complex)
    for i in range (NCOL):
        matrix[i][i] = complex(1.0,0.0)
    return matrix

def copy_fields(b):
    for j in range(NMAT):
        X_bak[j] = b[j]
    return X_bak

def rejected_go_back_old_fields(a):
    for j in range(NMAT):
        X[j] = a[j]
    return X

def refresh_mom():
    for j in range (NMAT):
        mom_X[j] = random_hermitian()
    return mom_X

def random_hermitian(): 
    tmp = np.zeros((NCOL, NCOL), dtype=complex)
 
    for i in range (NCOL):
        for j in range (i+1, NCOL): 

            r1, r2 = box_muller()
            tmp[i][i] = complex(0.0, 0.0)
            tmp[i][j] = complex(r1, r2)/math.sqrt(2)
            tmp[j][i] = complex(r1, -r2)/math.sqrt(2)
    return tmp

def hamil(X,mom_X):
    ham = potential(X)
    for j in range (NMAT):
        ham += 0.50 * np.trace(np.dot(mom_X[j],mom_X[j])).real
    return ham
         
def potential(X):
    s1, s2 = 0.0, 0.0
    for i in range (NMAT):
        s1 += 0.50 * np.trace(X[i] @ X[i])   
        s1 += (g/4.0)* np.trace(X[i] @ X[i] @ X[i] @ X[i])

    co = np.dot(X[0],X[1]) - np.dot(X[1],X[0])
    tr = np.trace(np.dot(co,co))
    s2 -= 0.50 * h *(tr.real)
    return ((s1+s2).real)*NCOL

def force(X):
    f_X[0] = X[0] + g*(X[0] @ X[0] @ X[0])
    f_X[0] -= h*dagger(comm(X[1], comm(X[0], X[1])))
    f_X[1] = X[1] + g*(X[1] @ X[1] @ X[1])
    f_X[1] -= h*dagger(comm(X[0], comm(X[1], X[0])))
    return (f_X)*NCOL

def leapfrog(X,dt):

    mom_X = refresh_mom()
    ham_init = hamil(X,mom_X)

    for j in range(NMAT):
        X[j] += mom_X[j] * dt * 0.50

    for j in range(NMAT):

        for i in range(1, nsteps):
            f_X = force(X)
            mom_X[j] -= f_X[j]*dt
            X[j] += mom_X[j]*dt

    f_X = force(X)

    for j in range(NMAT):

        mom_X[j] -= f_X[j] * dt
        X[j] += mom_X[j] * dt * 0.50

    ham_final = hamil(X,mom_X)
    
    return X, ham_init, ham_final

def update(X,acc_count):

    X_bak = copy_fields(X)
    X, start, end = leapfrog(X, dt)
    change = end - start
    expDS.append(np.exp(-1.0*change))
    if np.exp(-change) < random.uniform(0,1):
        X = rejected_go_back_old_fields(X_bak)
        print(("REJECT: deltaS = " "%8.7f " " startS = " "%8.7f" " endS = " "%8.7f" % (change, start, end)))
    else:
        print(("ACCEPT: deltaS = " "%8.7f " "startS = " "%8.7f" " endS = " "%8.7f" % (change, start, end)))
        acc_count += 1 

    w, v = LA.eigh(X[0])
    w1, v1 = LA.eigh(X[1])
    evals.append(w/NCOL)
    evals.append(w1/NCOL)

    tmp0 = np.trace(X[0] @ X[0]).real  
    tmp1 = np.trace(X[1] @ X[1]).real 
    tmp2 = (tmp0+tmp1)*0.5/NCOL
    trX2.append(tmp2) 

    tmp00 = np.trace(X[0] @ X[0] @ X[0] @ X[0]).real 
    tmp11 = np.trace(X[1] @ X[1] @ X[1] @ X[1]).real
    tmp22 = (tmp00+tmp11)*0.5/NCOL
    trX4.append(tmp22)

    if MDTU%GAP == 0:
        f3.write("%4.8f  \t %4.8f \n" %(tmp0/NCOL, tmp1/NCOL))
        f4.write("%4.8f  \t %4.8f \n" %(tmp00/NCOL, tmp11/NCOL))

    return X,acc_count 

# The main routine 

if __name__ == '__main__':
    
    
    if READIN == 0:
        for i in range (NMAT):
            for j in range (NCOL):
                for k in range (NCOL):
                    X[i][j][k] = complex(0.0,0.0)
                    
    if READIN == 1:
        print ("Reading old config.")
        name_f = "config_2MM_N{}.txt".format(NCOL)
        with open(name_f) as f2:
            A = np.loadtxt(f2).view(complex)
        f2.close()

        for i in range (NMAT):
            for j in range (NCOL):
                for k in range (NCOL):
                    X[i][j][k] = A[(NCOL*i)+j][k]

    f3 = open("t2.txt", "w")
    f4 = open("t4.txt", "w")

    acc_count = 0.0
    for MDTU in range (Niters_sim):
        X,acc_count = update(X,acc_count)

        if (MDTU+5)%10 == 0 and SAVE == 1:

            print ("Saving config.")
            name_f = "config_2MM_N{}.txt".format(NCOL)
            f1 = open(name_f, "w")
            for i in range (NMAT):
                np.savetxt(f1, X[i].view(float), delimiter= " ")  
            f1.close()

    f3.close()
    f4.close()
   
    if READIN == 0:
        expDS = expDS[cut:] 

    if acc_count/Niters_sim < 0.50:
        print("WARNING: Acceptance rate is below 50%")

    t2t4_plot = plt.figure(1) 
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(which='major', axis='y', linestyle='--')
    MDTU = np.linspace(0, int(Niters_sim/GAP), int(Niters_sim/GAP), endpoint=True)
    plt.xlabel('Time units')
    plot(MDTU, trX2, 'blue', label=r'Tr(X$^2$)')
    plot(MDTU, trX4, 'red', label=r'Tr(X$^4$)')
    plt.axhline(y=0.421783612, color='blue', linestyle='--')
    plt.axhline(y=0.333341358, color='red', linestyle='--')
    # Bootstrap results
    plt.legend(loc='best')
    evals = np.reshape(evals, (NCOL*Niters_sim*NMAT)) 
    hist_plot = plt.figure(2)
    plt.hist(evals, density=False, bins=50) 
    plt.ylabel(r'$\rho$')
    plt.xlabel('Eigenvalues (normalized by 1/N)')
    pp = PdfPages("2MM_allplots.pdf")
    pp.savefig(t2t4_plot, dpi = 300, transparent = True)
    pp.savefig(hist_plot, dpi = 300, transparent = True)
    pp.close()
    print("<Tr X^2 / NCOL>", np.mean(trX2), "+/-", (np.std(trX2)/np.sqrt(np.size(trX2) - 1.0)))
    print("<Tr X^4 / NCOL>", np.mean(trX4), "+/-", (np.std(trX4)/np.sqrt(np.size(trX4) - 1.0)))
    print("exp(-deltaH)", np.mean(expDS), "+/-", np.std(expDS)/np.sqrt(np.size(expDS) - 1.0))
    print ("COMPLETED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S")) 
    endTime = time.time() 
    print ("Running time:", round(endTime - startTime, 2),  "seconds")
    
