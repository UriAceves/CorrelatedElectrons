import  math
import  numpy  as np
class  Grid:
    # Constructor
    def  __init__(self , Z, rmin , rmax , dx):
        self.Z   = Z
        # Atomic  number
        xmin     = math.log(Z*rmin)
        xmax     = math.log(Z*rmax)
        self.dx = dx
        # Delta x
        self.x   = np.arange(xmin , xmax , dx)
        # Transformed  uniform  grid
        self.r   = np.exp(self.x)/Z
        # Logarithmic  grid
        self.V   = -Z/self.r
        # Potential  on the  grid
        self.N   = len(self.r)
        # Number  of grid  points

    # Print  out  grid  information
    def  __str__(self):
        info =   "==============================\n"
        info += "Grid  information :\n"
        info += "         Z: %18d\n"     % self.Z
        info += "     rmin: %18.10f\n" % self.r[0]
        info += "     rmax: %18.10f\n" % self.r[-1]
        info += "       dx: %18.10f\n" % self.dx
        info += "         N: %18d\n"     % self.N
        info += "==============================\n"
        return  info
