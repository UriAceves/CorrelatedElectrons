import  math
import  numpy  as np
import scipy.integrate as sp
class  Orbit:
    # Constructor
    def  __init__(self , G, n, l, occ=0, isSolved=True):
        self.G     = G        # Grid
        self.n     = n        # Principal  quantum  number
        self.l     = l        # Angular  momentum  quantum  number
        self.occ   = occ        # Number  of  occupied  electrons
        self.R = np.zeros(G.N)        # Radial  wave  function
        self.u = np.zeros(G.N)        # u = r*R
        self.E     = 0        # Orbital  energy
        self.node = 0        # Number  of  nodes
        self.M     = 0        # Matching  point  index
        self.dE    = 0        # Proposed  dE from  purturbation
        self.Vh = np.zeros(G.N)
        self.vreal = np.zeros(G.N)
        if  isSolved: self.solve() # Return  the  solved  orbital
        if isSolved: self.get_vh()

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
        (N, dx , r, V, l, u, Vh) = (
        self.G.N, self.G.dx, self.G.r, self.G.V, self.l, self.u, self.Vh)
        # Matching  point
        M = np.where(E<V)[0][0]  if E<V[-1] else N//2

        # Forward  initial  condition
        u[:2] = r[:2]**(l+1)/np.sqrt(r[:2])
        # Backward  initial  condition
        end = N-1
        if E<V[-1]:
            tail = np.exp(-math.sqrt(abs(2*E))*r)/np.sqrt(r)
            end   = np.where(tail >0) [0][ -1]
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
        # Propose  dE by  purturbation
        kku = 6/5*(   2/dx**2*u[M]-(1/dx**2+kk[M+1]/12)*u[M+1] -(1/dx**2+kk[M -1]/12)*u[M-1] )
        dE   = 0.5*(kku -kk[M]*u[M])*u[M]*dx
        # Count  the  number  of nodes
        node = np.sum(u[:-1]*u[1:] <0)
        # Collect  results
        (self.R, self.u, self.E, self.node , self.M, self.dE) = (u/np.sqrt(r), u*np.sqrt(r), E, node , M, dE)

    # Shoot  repeatly  within (Em, Ep) until  eigen -energy
    def  shootbracket(self , Em , Ep):
        # Collect  fields
        (r, N) = (self.G.r, self.G.N)
        # Error  message
        errmsg = "Orbit: not  found  in (%.2f, %.2f)\n" % (Em , Ep)
        # Determine  the  boundary  dE
        self.shoot(Em); dEm = self.dE
        self.shoot(Ep); dEp = self.dE
        if dEm*dEp >0:
            print(errmsg); return  False
        itMax = 100
        for it in  range(itMax):
            E = self.E+self.dE
            # Prevent  going  beyond  bracket
            if E<=Em or E>=Ep:
                E = 0.5*(Em+Ep)
            # Determine  the dE of  current  energy
            self.shoot(E); dE = self.dE
            # Check  dE
            if abs(dE)<1e-8:
                return  True
            # Converged
            else:
                (Em , Ep) = (E, Ep) if dE*dEm >0 else (Em, E)
        # Cannot  converge  to the  tolerance
        print(errmsg); return  False

    # Solve  the eigen -state
    def  solve(self):
        # Collect  fields
        (Z, n, l) = (self.G.Z, self.n, self.l)
        # Big  bracket
        (Em , Ep) = ( -1.1*0.5*Z**2,  0.1*0.5*Z**2)
        # Error  message
        errmsg = "Orbit: n=%d l=%d not  found  in (%.2f, %.2f)\n" % (n, l, Em, Ep)
        # Check  lower  and  upper  bound
        self.shoot(Em); (ndm , dEm) = (self.node , self.dE)
        self.shoot(Ep); (ndp , dEp) = (self.node , self.dE)
        if ndm >n-l-1 or ndp <n-l-1 or (ndm==n-l-1 and dEm <0) or (ndp==n-l-1 and dEp >0):
            print(errmsg); return  False
        itMax = 100
        mPass = pPass = False
        # Find  lower  bracket
        (em , ep) = (Em , Ep)
        for it in  range(itMax):
            E = 0.5*(em+ep)
            self.shoot(E); (nd , dE) = (self.node , self.dE)
            if nd==n-l-1 and dE >0:
                (Em , mPass) = (E, True); break
            else:
                (em , ep) = (E, ep) if nd <n-l-1 else (em, E)
        # Find  upper  bracket
        (em , ep) = (Em , Ep)
        for it in  range(itMax):
            E = 0.5*(em+ep)
            self.shoot(E); (nd , dE) = (self.node , self.dE)
            if nd==n-l-1 and dE <0:
                (Ep , pPass) = (E, True); break
            else:
                (em , ep) = (em , E) if nd >n-l-1 else (E, ep)
        if  mPass  and  pPass: # Bracket  found
            return  self.shootbracket(Em, Ep)
        else:
            print(errmsg); return  False

    def get_vh(self):
        (N, dx ,r, u, Vh) = (self.G.N, self.G.dx, self.G.r, self.u, self.Vh)
        #Charge as function of r
        Q = np.zeros(N)
        #Uniform charge
        Q = sp.cumtrapz(r**3*(3./r[-1]**3) ,dx=dx,initial=0)
        # Vh = [np.trapz(Q[i:]/r[i:], dx=dx) + 1./r[-1] for i in range(0,N)]
        # vreal = [1./(2*r[-1]) + -r[i]**2/(2*r[-1]**3) + 1/r[-1] for i in range(0,N)]
        #Q = sp.cumtrapz(u**2*r ,dx=dx,initial=0)
        Vh = [np.trapz(Q[i:]/r[i:], dx=dx) + 1./r[-1] for i in range(0,N)]
        vreal = [1./(2*r[-1]) + -r[i]**2/(2*r[-1]**3) + 1/r[-1] for i in range(0,N)]
        #Hartree potential
        self.vreal = vreal
        self.Vh = Vh
        #return Vh
