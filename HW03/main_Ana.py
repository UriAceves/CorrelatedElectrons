from grid import Grid
from orbit import Orbit
import matplotlib.pyplot as plt
import numpy as np


Z = 1
g = Grid(Z, 0.0001/Z, 110., 0.001)
n = 1

epsilon = 0.0001

colors = ["black","red","blue","green","purple","orange"]
styles = ["-","-.","-","--"]

if n<=6:
	factor = 0.3
else:
	factor = 0.1 

for l in range(0,n):
	orb = Orbit(g, n, l)
	E = -0.45
	orb.shoot(E)
	if n-l-1==0:
	  while abs(orb.delta_E)>epsilon:
		 E += 1e-4*orb.delta_E
		 orb.shoot(E)
		 print "----->"+str(E)
	else:
	  while abs(orb.delta_E)>epsilon or orb.node != n-l-1:
		 E += 1e-4*orb.delta_E
		 if abs(orb.delta_E)<epsilon:
		   E -= factor*E
		 print "--------->"+str(E)
		 orb.shoot(E)
	#print "I have done l="+str(l)
	plt.plot(g.r, orb.u, linewidth = 2.0, color=colors[l%6], linestyle = styles[l%4], label = "$u_{"+str(n)+","+str(l)+"}$")  

#print "El numero de nodos es "+str(orb.node)+"---->\tn="+str(n-l-1)+"\t"+str(E)
#print orb.tp

ax = [0 for i in g.r]

plt.plot(g.r, ax, color="black")
#plt.plot(g.r, energy, color="red")
plt.xlabel("$r$")
plt.ylabel("u")
plt.ylim(-0.2,0.8)
plt.xlim(0,20)
plt.legend()

plt.show()
