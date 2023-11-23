import numpy as np

from .abstract_dynamic_system import AbstractDynamicSystem

g, l = 9.8, np.sqrt(0.5**2 + 1**2)

def F(X, t, LAMBDA, w0):
	return np.array([X[1], -2*LAMBDA*X[1]-w0**2*X[0]],np.float64)

## Dummy dynamic system just to test
class RessortTorsionDynamicSystem(AbstractDynamicSystem):

    def __init__(self, mesh):
        ## Constructor
        # @param self
        # @param mesh  
        # @param theta0  
        # @param lambda  
        # @param w0  
        super().__init__()
        self.mesh = mesh

        # Animations parameters
        self.pulsation_propre = 10
        self.coef_amortissement = 1
        self.it = 60.
        self.delta = 0.01
        self.period = 120.
        self.theta0 = np.arctan(self.mesh.positions[1]/self.mesh.positions[2])
        self.X = np.array([self.theta0,0],np.float64) # [theta, theta']
        

    def step(self, method="explicit"):

        if method=="explicit":
            self.X += self.delta * F(self.X, self.it, self.coef_amortissement, self.pulsation_propre)
        elif method=="implicit":
            pass
        theta = self.X[0]
        self.mesh.positions[2], self.mesh.positions[3] = l*np.sin(theta), l*(1-np.cos(theta))
        self.it += self.delta