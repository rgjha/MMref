import numpy  as np
import matplotlib.pyplot as plt

x1,x2,x3,x4,x5 = np.loadtxt("t2", delimiter='\t', unpack=True)
y1,y2,y3,y4,y5 = np.loadtxt("t4", delimiter='\t', unpack=True)
z = np.loadtxt("action_N300_D5.txt")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(which='major', axis='y', linestyle='--')
plt.xlabel('Time units')
plt.ylabel(r'$\langle S/(N^2-1) \rangle $')
plt.plot(z)
'''
plt.plot(x1)
# ,label=r'$\langle Tr(X^2)\rangle$'
plt.plot(x2)
plt.plot(x3)
plt.plot(x4)
plt.plot(x5)
plt.plot(y1)
plt.plot(y2)
plt.plot(y3)
plt.plot(y4)
plt.plot(y5)
'''
#plt.plot(x3,'green',label=r'$t_{4}$')
#plt.plot(x4,'orange',label=r'$t_{4, II}$')
#plt.plot(y,'blue', label=r'$t_{4}$, closed')
#plt.plot(z,'teal', label=r'$t_{2}$, open')
#plt.plot(w,'green', label=r'$t_{4}$, open')
plt.legend(loc='best')
plt.savefig('act_D5_YM.pdf')
plt.show()