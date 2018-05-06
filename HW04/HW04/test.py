from  grid  import  Grid
from  orbit  import  Orbit
import  matplotlib.pyplot  as plt
Z = 1
G = Grid(Z, 0.0001/Z, 50, 0.002)
orb = Orbit(G,3, 0)
# That’s it!
# Print  resulting E
#print(orb.E)
#print(orb.Vh)
# Plot
#plt.plot(G.r, orb.u, linewidth =2)
#plt.title(r"Numerical solution $u_{30}$")
plt.plot(G.r, orb.vreal, "y-", linewidth =4, label="Analytical")
plt.plot(G.r, orb.Vh, "k--",linewidth =1, label="Numerical")
plt.legend(loc="best")
plt.savefig("num30.png")
plt.show()
