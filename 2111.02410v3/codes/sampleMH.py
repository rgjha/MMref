import random
import math
import time
import datetime
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


startTime = time.time()
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S"))

if len(sys.argv) < 2:
    print("Usage:python",str(sys.argv[0]),"exp(-x) [option 1] or exp(-x**2) [option 2]")
    sys.exit(1)

option = int(sys.argv[1])

if option not in [1,2]:
    print ("Wrong input. ABORT")
    sys.exit(1) 

niter = int(1e+6) 
naccept = 0
x = []

def target(input):

    if option == 1:
        if input < 0: 
            return 0
        else:
            return np.exp(-input)  

    if option == 2: 
        return np.exp(-input**2)  


for iter in range(1,niter):

    if iter == 1: 
        current_x = 0.0

    choose = np.random.normal(0,1)
    proposed_x = current_x + choose
    A = min(1, target(proposed_x)/target(current_x)) 

    if random.uniform(0,1) < A:
        x.append(proposed_x)
        naccept += 1 
        current_x = proposed_x
    else:
        x.append(current_x)

print ("Accepted # ", naccept)


if option == 1:
    n, bins, patches = plt.hist(x, 30, facecolor='green', alpha=0.55, density=True, range=[0, 4])
    plt.xlabel('x')
    plt.ylabel('Probability')
    #plt.title(r'$\mathrm{Histogram\ of\ data sampled by Metropolis:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
    plt.grid(True)
    plt.savefig('exp(-x)metro.pdf')

if option == 2:
    # If we were fitting a gaussian 
    (mu, sigma) = norm.fit(x)
    # Histogram, denisty=True gives probability not frequency 
    n, bins, patches = plt.hist(x, 30, facecolor='green', alpha=0.55, density=True, range=[-4, 4])
    # Best fit line (if needed)
    y = norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('Probability')
    #plt.title(r'$\mathrm{Histogram\ of\ data sampled by Metropolis:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
    plt.grid(True)
    plt.savefig('exp(-x2)metro.pdf')
