from  grid  import  Grid
from  orbit  import  Orbit
import  matplotlib.pyplot  as plt
# Set up grid
Z = 1
G = Grid(Z, 0.0001/Z, 80, 0.002)
# Set up  orbit
n = 2
l = 0
expected_nodes = n-l-1
orb = Orbit(G, n, l)

# Shoot  with a guess E
E =  -0.8
orb.shoot(E)
eps = 0.0001
if n-l-1 ==0:
    while abs(orb.delt_E) > eps:
        E += orb.delt_E*1e-4
        orb.shoot(E)
        print(E)
else:
    while abs(orb.delt_E) > eps or orb.node != n-l-1:
        E += orb.delt_E*1e-4
        if abs(orb.delt_E) < eps:
            E -= .3*E
        orb.shoot(E)
        print(E)

print("Number of nodes ",orb.node)
#print(orb.tp)
# Plot
plt.plot(G.r, orb.u, linewidth =2.0)
plt.xlim(0, 80)
plt.ylim(-2, 2)
plt.xlabel(r"$\ (\mathrm{a_0})$", fontsize =16)
plt.ylabel(r"$u_{nl}$", fontsize =16)
print(orb)
plt.show()
