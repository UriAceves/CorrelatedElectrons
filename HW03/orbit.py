import  math
import  numpy  as np
class  Orbit:
    # Constructor
    def  __init__(self , G, n, l, occ=0):
        self.G     = G  # Grid
        self.n     = n  # Principal  quantum  number
        self.l     = l  # Angular  momentum  quantum  number
        self.occ   = occ    # Number  of  occupied  electrons
        self.R = np.zeros(G.N)  # Radial  wave  function
        self.u = np.zeros(G.N)  # u = r*R
        self.E     = 0  # Orbital  energy
        self.node = 0   # Number  of  nodes
        self.M     = 0  # Matching  point  index
        self.tp = np.zeros(2) #Classical turning points
        self.delt_E = 0. #Energy change

    # Print  out  orbital  information
    def  __str__(self):
        info =  "==============================\n"
        info += "Orbit  information :\n"
        info += "         E: %18.10f\n" % self.E
        info += "     node: %18d\n"     % self.node
        info += "         n: %18d\n"     % self.n
        info += "         l: %18d\n"     % self.l
        info += "      occ: %18d\n"     % self.occ
        info += "         M: %18d\n"     % self.M
        info += "==============================\n"
        return  info

    # Shoot a wave  function  with  given E
    def  shoot(self , E):
        # Collect  fields
        (N, dx , r, V, l, u, tp) = (
        self.G.N, self.G.dx, self.G.r, self.G.V, self.l, self.u, self.tp)

        #Calculating the classical turning points
        tp[0] = (-self.G.Z + np.sqrt(self.G.Z**2+2*E*l*(l+1)))/(2*E)
        tp[1] = (-self.G.Z - np.sqrt(self.G.Z**2+2*E*l*(l+1)))/(2*E)

        #print("The turning points are ", turning_points[0], turning_points[1])

        # Matching  point
        M = np.where(E<V)[0][0]  if E<V[-1] else N//2

        # Forward  initial  condition
        u[:2] = r[:2]**(l+1)/np.sqrt(r[:2])
        # Backward  initial  condition
        end = N-1
        if E<V[-1]:
            tail = np.exp(-math.sqrt(abs(2*E))*r)/np.sqrt(r)
            end   = np.where(tail >0) [0][-1]
            u[end -1:] = tail[end -1:]
        else:
            # E too high , then  particle  in a box
            u[-2:] = [1e-6, 0.0]

        # Numerov  method
        kk = 2*(r**2*(E-V) - 0.5*(l+0.5) **2)
        A   = 2 - dx*dx*5/6* kk
        B   = 1 + dx*dx/12*kk
        # Forward  integration
        for i in range(2, M+1):
            u[i] = (A[i-1]*u[i-1] - B[i-2]*u[i-2]) / B[i]
        FM = u[M]

        # Backward  integration
        for i in range(end -2, M-1,  -1):
            u[i] = (A[i+1]*u[i+1] - B[i+2]*u[i+2]) / B[i]
        BM = u[M]

        # Connect  forward  and  backward  parts
        u[M:] *= FM/BM
        # Normalization
        u /= math.sqrt(np.trapz(r**2*u**2, dx=dx))

        #Calculating Delta k to get the next energy change
        k_num = ((u[M+1]*B[M+1] + u[M-1]*B[M-1])/u[M] - 2)*(-6./(5*dx*dx))
        print(k_num -kk[M])
        E_num = (k_num/2 + r[M]**2*V[M] + 0.5*(l+0.5)**2)/r[M]**2
        self.delt_E = E_num - E

        #Contributions from the missing regions
        integral = 0.5*((u[0]*r[0])**2*self.G.x[0])

        if tail[-1] != 0:
            rad = 1.1*r[-1]
            while np.exp(-math.sqrt(abs(2*E))*rad)/np.sqrt(rad) > 0.:
                rad = 1.1*rad
            integral += 0.5*((u[-1]*r[-1])**2)*(np.log(self.G.Z*rad)-
                                                    np.log(self.G.Z*r[-1]))

        print("The area of the missing regions is approximately ", integral)

        # Count  the  number  of nodes
        node = np.sum(u[:-1]*u[1:] <0)
        # Collect  results
        (self.R, self.u, self.E, self.node , self.M, self.tp) = (u/np.sqrt(r), u*np.sqrt(r), E, node , M, tp)
