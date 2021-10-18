import numpy  as np
import matplotlib.pyplot as plt

#x1 = np.loadtxt("t2_N300_closed", delimiter='\t', unpack=True)
#x2 = np.loadtxt("t4_N300_closed", delimiter='\t', unpack=True)

x1 = np.loadtxt("t2_g1")
x2 = np.loadtxt("t4_g1")
x3 = np.loadtxt("t2_g2")
x4 = np.loadtxt("t4_g2")


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(which='major', axis='y', linestyle='--')
plt.xlabel('Time units')
#plt.ylabel(r'$t_{2}, open$')
plt.plot(x1,'red',label=r'$t_{2}, g=1$')
plt.plot(x2,'blue',label=r'$t_{4}, g=1$')
plt.plot(x3,'green',label=r'$t_{2}, g=2$')
plt.plot(x4,'orange',label=r'$t_{4}, g=2$')
#plt.plot(x3,'green',label=r'$t_{4}$')
#plt.plot(x4,'orange',label=r'$t_{4, II}$')
#plt.plot(y,'blue', label=r'$t_{4}$, closed')
#plt.plot(z,'teal', label=r'$t_{2}$, open')
#plt.plot(w,'green', label=r'$t_{4}$, open')
plt.legend(loc='best')
plt.savefig('4MM_g1_g2.pdf')
plt.show()